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

use crate::cli::run_output::size_of_ground;
use crate::fasta::ContigList;
use crate::reference_info::ReferenceInfo;
use crate::region_typing::segment_criteria::{MinCopies, SsrSegmentCriteria};
use crate::region_typing::{GenomeRegions, RegionKind, TypedRegion, TypedRegionConfig};
use crate::regions::{BedError, ContigBounds, ContigTooLong};
use crate::repeat_catalog::{
    CriteriaRefusal, ReadScope, RepeatCatalog, RepeatCatalogError, StrRepeatCriteria,
    sibling_catalog_path,
};
use crate::run::report::{describe, name_of};
use crate::run::{
    RepeatGround, RepeatGroundSetting, RepeatSorting, RunError, Segmentation, SegmentationInputs,
};
use crate::tandem_repeat::{PeriodRange, PeriodRangeError};
use crate::types::{Bp, ContigId};

/// **What this run counts as a repeat** — the six flags that say so, as the three subcommands
/// that type spell them.
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
    /// An STR bundle whose hull is longer than this many bases is a satellite.
    pub max_str_bundle_len: u64,
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
    /// What this run counts as a repeat. **Owned**, because it is six small values a
    /// subcommand builds from its flags rather than something it already holds.
    pub routing: RepeatRouting,
    /// Which path each kind of repeat ground goes to (`--call-satellites`,
    /// `--generic-path-only`).
    pub setting: RepeatGroundSetting,
}

impl GroundRequest<'_> {
    /// The catalog this request names — the flag's path, or the one beside the reference.
    #[must_use]
    pub fn catalog_path(&self) -> PathBuf {
        catalog_path_for(self.catalog, self.reference)
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
        /// The file whose contig table it is: the reference, or a psp.
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
        "no repeat catalog at {}; build one with `pvc repeat-catalog --reference {}`",
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

    /// Both `--call-satellites` and `--generic-path-only` were given (`command_line.md`
    /// §4.4.1). The second already sends satellites to the generic path, so the pair means
    /// nothing of its own, and the refusal names the one to keep for each intent.
    #[error(
        "--generic-path-only already sends satellites to the SNP/indel model, so it is never \
         given with --call-satellites: keep --generic-path-only to call every repeat that way, \
         or --call-satellites to keep the repeat model for repeats up to --max-str-len"
    )]
    BothSettingsForRepeatGround,

    /// `call-from-psps --regions` asked for ground the psps hold no records for.
    ///
    /// **Refused before anything is called, and before the reference is read**: a psp holds
    /// records only for the ground its walk covered, so a run over more would report the rest
    /// as if it had been looked at and found empty. The message lists each piece outside, what
    /// the psps hold on that piece's contig, and the psps' whole ground.
    #[error("{}", message_for_regions_outside_the_psps(.bed, .pieces, .unlisted, .psps_ground))]
    RegionsOutsideThePsps {
        /// The BED given to `--regions`.
        bed: PathBuf,
        /// The first pieces outside, at most [`PIECES_LISTED`], in genomic order.
        pieces: Vec<RegionOutsideThePsps>,
        /// How many more pieces there are past those listed.
        unlisted: usize,
        /// The ground the psps were written over, named and sized for one line.
        psps_ground: String,
    },

    /// A catalog was named for psps written with `--generic-path-only`, which read none and
    /// sorted no repeats. Reading one now could only disagree with them.
    #[error(
        "--catalog {} was given, but the psps were written with --generic-path-only, which reads \
         no catalog; leave --catalog out",
        path.display()
    )]
    CatalogUnderGenericPathOnly {
        /// The catalog named.
        path: PathBuf,
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
    let bounds = addressable_bounds(request.reference, contigs)?;
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

/// How many pieces outside the psps a refusal lists before it says how many more there are.
pub const PIECES_LISTED: usize = 10;

/// One piece of `--regions` that lies outside the psps' ground, with what the psps hold on its
/// contig — both already written as a person reads them, 1-based and inclusive.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RegionOutsideThePsps {
    /// The piece, as `chr2:900001-1000000`.
    pub piece: String,
    /// What the psps hold on the piece's contig: `the psps cover nothing on chr3`, or
    /// `on chr2 the psps hold chr2:1-900000`.
    pub what_the_psps_hold_there: String,
}

/// **The regions `call-from-psps --regions` asks for, read from the BED and checked to lie inside
/// the ground the psps were written over.**
///
/// The BED is resolved against the psps' own contig table, so this needs nothing from the
/// reference and a bad BED is refused in seconds. It is read the way every `--regions` is read —
/// sorted, with overlapping and touching lines merged, and a line running past its contig's end
/// cut at that end — and then every base of it must lie in `psps_ground`. `first_psp` names the
/// file the contig table came from, for the one refusal about the table itself.
///
/// # Errors
///
/// [`GroundError::Bed`] for a BED that will not read — a contig the psps' table does not have, a
/// malformed line, or no interval at all; [`GroundError::RegionsOutsideThePsps`] for one that asks
/// for any base the psps were not written over; and [`GroundError::ContigTooLong`] for a table a
/// run cannot address.
pub fn regions_inside_the_psps(
    bed: &Path,
    first_psp: &Path,
    psps_contigs: &ContigList,
    psps_ground: &GenomeRegions,
) -> Result<GenomeRegions, GroundError> {
    let bounds = addressable_bounds(first_psp, psps_contigs)?;
    let requested =
        GenomeRegions::from_bed_path(bed, &bounds).map_err(|source| GroundError::Bed {
            path: bed.to_path_buf(),
            source,
        })?;
    let outside = requested.parts_outside(psps_ground);
    if outside.is_empty() {
        return Ok(requested);
    }
    let pieces = outside
        .iter()
        .take(PIECES_LISTED)
        .map(|piece| RegionOutsideThePsps {
            piece: name_of(*piece, psps_contigs),
            what_the_psps_hold_there: what_the_psps_hold_on(
                piece.contig,
                psps_contigs,
                psps_ground,
            ),
        })
        .collect();
    let bases: u64 = psps_ground.iter().map(|region| region.len()).sum();
    Err(GroundError::RegionsOutsideThePsps {
        bed: bed.to_path_buf(),
        pieces,
        unlisted: outside.len().saturating_sub(PIECES_LISTED),
        psps_ground: format!(
            "{} ({})",
            describe(psps_ground, psps_contigs),
            size_of_ground(bases)
        ),
    })
}

/// What the psps hold on one contig, as the second column of the refusal writes it.
fn what_the_psps_hold_on(contig: ContigId, contigs: &ContigList, ground: &GenomeRegions) -> String {
    let contig_name = contigs
        .entries
        .get(contig.get() as usize)
        .map_or("this contig", |entry| entry.name.as_str());
    let on_this_contig = ground.on_contig(contig);
    if on_this_contig.is_empty() {
        return format!("the psps cover nothing on {contig_name}");
    }
    let on_this_contig = describe(&on_this_contig, contigs);
    format!("on {contig_name} the psps cover {on_this_contig}")
}

/// The refusal's text: what was asked, each piece outside beside what the psps hold on its contig,
/// the psps' whole ground, and what to do.
fn message_for_regions_outside_the_psps(
    bed: &Path,
    pieces: &[RegionOutsideThePsps],
    unlisted: &usize,
    psps_ground: &str,
) -> String {
    let width = pieces
        .iter()
        .map(|piece| piece.piece.len())
        .max()
        .unwrap_or(0);
    let mut lines = vec![format!(
        "--regions {} asks for regions the psps do not cover, so nothing was called:",
        bed.display()
    )];
    lines.extend(pieces.iter().map(|piece| {
        format!(
            "  {:<width$}   {}",
            piece.piece, piece.what_the_psps_hold_there
        )
    }));
    if *unlisted > 0 {
        lines.push(format!("  … and {unlisted} more"));
    }
    lines.push(format!(
        "The psps cover {psps_ground}. Positions are 1-based, both ends included."
    ));
    lines.push(
        "Remove these regions from the BED, or write the psps again over regions that include them."
            .to_owned(),
    );
    lines.join("\n")
}

/// The contig table as the BED reader takes it, **refusing a contig longer than a run can
/// address** rather than resolving against a wrong length.
///
/// `ContigBounds` carries a `u32`, and a contig past that would have its ground resolved against a
/// wrong length with nothing to notice — the rule `typed_regions.rs`'s `ContigTooLong` records, and
/// the reason the one precedent in this tree that casts is a test. `source` is the file the table
/// came from, named by the refusal: the reference, or a psp.
fn addressable_bounds<'a>(
    source: &Path,
    contigs: &'a ContigList,
) -> Result<Vec<ContigBounds<'a>>, GroundError> {
    // **The narrowing is refused, not taken** (`ContigBounds::narrowed` says why).
    ContigBounds::of_contigs(
        contigs
            .entries
            .iter()
            .map(|entry| (entry.name.as_str(), entry.length)),
    )
    .map_err(
        |ContigTooLong { name, length }| GroundError::ContigTooLong {
            reference: source.to_path_buf(),
            name,
            length,
            limit: u64::from(u32::MAX),
        },
    )
}

/// **The catalog a command reads**: the path it was given, or the file `repeat-catalog` writes
/// beside the reference.
///
/// **One rule, in one place.** It is part of what a user is promised — every subcommand's
/// `--catalog` says *"defaults to `<reference>.repeats.parquet`"* — and it is spelled at four
/// call sites, with a fifth arriving with `regenerate-census`. A copy that drifts is a command
/// that silently reads a different file from its siblings.
#[must_use]
pub fn catalog_path_for(catalog: Option<&Path>, reference: &Path) -> PathBuf {
    catalog.map_or_else(|| sibling_catalog_path(reference), Path::to_path_buf)
}

/// **What this run counts as a repeat**, built from the six flags that say so.
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
        max_str_bundle_len: Bp(routing.max_str_bundle_len),
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
/// from the catalog, **with what counts as a repeat taken from this run's six flags**.
///
/// A command whose criteria come from somewhere else — a psp header, which carries the finished
/// value the walk routed under — calls [`segments_cut_with`] instead, and everything past the
/// flags is shared.
///
/// **Under `--generic-path-only` no catalog is opened**, and the segments are
/// [`segments_without_a_catalog`]'s.
///
/// # Errors
///
/// [`GroundError::MissingCatalog`] first, then [`GroundError::PeriodRange`] when the two period
/// flags do not make a range, and then everything [`segments_cut_with`] answers with.
pub fn segments_over(
    request: &GroundRequest<'_>,
    analysed: &GenomeRegions,
    with_checksums: &ReferenceInfo,
) -> Result<Segmentation, GroundError> {
    if request.setting == RepeatGroundSetting::GenericPathOnly {
        return Ok(segments_without_a_catalog(analysed));
    }
    let path = request.catalog_path();
    // **Before the flags are converted**, so that a run whose catalog is not there is told that
    // and not which of its period flags it typed backwards. [`segments_cut_with`] checks the
    // same thing again, for the callers that enter there.
    refuse_a_catalog_that_is_not_there(&path, request.reference)?;
    let criteria = routing_criteria(&request.routing)?;
    segments_cut_with(
        &path,
        request.reference,
        request.setting,
        &criteria,
        CriteriaSource::TheCommandLineFlags,
        analysed,
        with_checksums,
    )
}

/// The setting for repeat ground the two options choose, refusing both at once
/// ([`GroundError::BothSettingsForRepeatGround`]). **Asked before anything is read**, so a run
/// given both stops at its command line.
///
/// # Errors
///
/// [`GroundError::BothSettingsForRepeatGround`] when both are given.
pub fn setting_of_the_options(
    call_satellites: bool,
    generic_path_only: bool,
) -> Result<RepeatGroundSetting, GroundError> {
    if call_satellites && generic_path_only {
        return Err(GroundError::BothSettingsForRepeatGround);
    }
    Ok(RepeatGroundSetting::from_options(
        call_satellites,
        generic_path_only,
    ))
}

/// **The segments of a run that sorts no repeats** (`--generic-path-only`): one generic region
/// per analysed span, and no catalog opened.
///
/// Every repeat is ordinary sequence to such a run, so there is nothing to cut the analysed
/// ground at; and a generic region must be maximal, since an indel across a region's edge is
/// called by neither side (`region_typing::fill_generic_gaps`). The spans are the BED's as
/// [`GenomeRegions`] normalised them, as generic ground is clipped to them in every run.
#[must_use]
pub fn segments_without_a_catalog(analysed: &GenomeRegions) -> Segmentation {
    let segments = analysed
        .iter()
        .map(|region| TypedRegion {
            region,
            kind: RegionKind::Generic,
        })
        .collect();
    Segmentation::of_segments(segments, analysed.clone(), RepeatGround::GenericPathOnly)
}

/// The setting, catalog and criteria a cohort of psps was walked with, which a command reading
/// them sorts its own ground with again — `None` for psps written with `--generic-path-only`,
/// which sorted no repeats and read no catalog ([`segments_without_a_catalog`]).
#[must_use]
pub fn sorting_of_the_psps(
    inputs: &SegmentationInputs,
) -> Option<(RepeatGroundSetting, &RepeatSorting)> {
    inputs
        .repeat_ground
        .sorting()
        .map(|sorting| (inputs.repeat_ground.setting(), sorting))
}

/// Refuse a `--catalog` named for psps written with `--generic-path-only`
/// ([`GroundError::CatalogUnderGenericPathOnly`]). A catalog left unnamed is not looked for.
///
/// # Errors
///
/// As above.
pub fn refuse_a_catalog_the_psps_never_read(catalog: Option<&Path>) -> Result<(), GroundError> {
    match catalog {
        Some(path) => Err(GroundError::CatalogUnderGenericPathOnly {
            path: path.to_path_buf(),
        }),
        None => Ok(()),
    }
}

/// **Where a run's repeat criteria came from** — which is what a catalog that cannot serve them
/// has to know before it can name something the reader is able to move.
///
/// A catalog is built below every calling floor so that a reader filters rather than re-scans
/// (`repeat_catalog.md` §4.1); asking it for tracts *below* what it holds is refused, because
/// those rows were never written. What to do about that refusal depends entirely on where the
/// numbers came from: a person who typed `--min-copies 3,3,3,3,3,3` should be told that flag,
/// and a person whose numbers came out of a psp header has no flag to move — their catalog is
/// not the one their psps were walked against, and that is what they have to fix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CriteriaSource {
    /// The repeat-routing flags of this command line ([`RepeatRouting`]).
    TheCommandLineFlags,
    /// The psp headers of the cohort being read, where the walk recorded them
    /// (`psp_census_pair.md` §6).
    ThePspHeaders,
}

/// The same segments, cut with criteria the caller already holds.
///
/// **What a psp-driven command enters at.** The six flags are one way to say what this run
/// counts as a repeat; a psp header is the other, and it carries the finished
/// [`StrRepeatCriteria`] with no way back to the flags it was built from
/// (`psp_census_pair.md` §6). Everything past that — refusing a catalog that is not there,
/// checking the one that is against the reference, cutting the ground and recording what it was
/// cut with — is the same work, and this is it.
///
/// `catalog_path` is the file itself, already resolved; `reference` is named only so a refusal
/// can say which assembly the missing catalog belongs to; `criteria_from` decides what a catalog
/// too coarse for these criteria tells the reader to change.
///
/// # Errors
///
/// [`GroundError::MissingCatalog`] when there is no catalog to read,
/// [`GroundError::Catalog`] when it will not read, was built on another reference, or cannot
/// serve criteria that came from a psp header,
/// [`GroundError::RoutingBelowCatalog`] when the flags of this command line asked for repeats it
/// does not hold, and [`GroundError::Segmentation`] when its segment stream fails part-way.
pub fn segments_cut_with(
    catalog_path: &Path,
    reference: &Path,
    setting: RepeatGroundSetting,
    criteria: &StrRepeatCriteria,
    criteria_from: CriteriaSource,
    analysed: &GenomeRegions,
    with_checksums: &ReferenceInfo,
) -> Result<Segmentation, GroundError> {
    let catalog = open_catalog(catalog_path, reference, with_checksums)?;
    segments_cut_from(
        &catalog,
        catalog_path,
        setting,
        criteria,
        criteria_from,
        analysed,
    )
}

/// **The catalog at `catalog_path`, opened and checked against the reference** — the first half
/// of [`segments_cut_with`], for a command that has something to ask the catalog's header before
/// any segment is cut.
///
/// `estimate-parameters` and `regenerate-census` are those commands: each compares the header with
/// the catalog its psps were walked with, and a catalog that is not theirs would otherwise be met
/// while cutting, as one that cannot serve the criteria — a message about the file's contents where
/// the fault is which file was named.
///
/// # Errors
///
/// [`GroundError::MissingCatalog`] when there is no catalog to read, and [`GroundError::Catalog`]
/// when it will not read, was written by another build of this program, or was built on another
/// reference.
pub fn open_catalog(
    catalog_path: &Path,
    reference: &Path,
    with_checksums: &ReferenceInfo,
) -> Result<RepeatCatalog, GroundError> {
    refuse_a_catalog_that_is_not_there(catalog_path, reference)?;
    RepeatCatalog::open_checking_against_reference(catalog_path, with_checksums).map_err(|source| {
        GroundError::Catalog {
            path: catalog_path.to_path_buf(),
            source,
        }
    })
}

/// **The segments cut from a catalog already open** — the second half of [`segments_cut_with`].
///
/// `catalog_path` is the file `catalog` was opened from, recorded in the segmentation and named by
/// a refusal; an open catalog does not hand its path out.
///
/// **`CriteriaSource::ThePspHeaders` is now defensive for both psp-taking commands**, each of
/// which proves the catalog's header is the one its psps were walked with before it cuts anything:
/// the criteria it cuts with are those psps' own, and a catalog that cut them once cuts them
/// again. What is left for that arm is a file whose header matches and whose rows do not.
///
/// # Errors
///
/// As [`segments_cut_with`], less [`GroundError::MissingCatalog`]. [`GroundError::Catalog`] still
/// arrives, for a catalog too coarse for criteria that came from a psp header, and
/// [`GroundError::RoutingBelowCatalog`] for the same asked by this command line's own flags.
///
/// # Panics
///
/// Given [`RepeatGroundSetting::GenericPathOnly`], which reads no catalog: such a run's segments
/// are [`segments_without_a_catalog`]'s, and every caller routes it there.
pub fn segments_cut_from(
    catalog: &RepeatCatalog,
    catalog_path: &Path,
    setting: RepeatGroundSetting,
    criteria: &StrRepeatCriteria,
    criteria_from: CriteriaSource,
    analysed: &GenomeRegions,
) -> Result<Segmentation, GroundError> {
    let spans: Vec<_> = analysed.iter().collect();
    let segments = catalog
        .genome_segments(criteria, ReadScope::Regions(&spans))
        .map_err(|source| match criteria_from {
            CriteriaSource::TheCommandLineFlags => {
                catalog_error_naming_the_flag(source, catalog_path)
            }
            // **No flag to name**: these criteria are what the walk recorded, and this command
            // line cannot change them. What the reader can change is which catalog they pointed
            // it at, which the general refusal already names.
            CriteriaSource::ThePspHeaders => GroundError::Catalog {
                path: catalog_path.to_path_buf(),
                source,
            },
        })?;
    let segments =
        segments
            .collect::<Result<Vec<_>, _>>()
            .map_err(|source| GroundError::Segmentation {
                source: RunError::Catalog {
                    path: catalog_path.to_path_buf(),
                    source,
                },
            })?;
    let sorting = RepeatSorting {
        catalog: catalog.header().clone(),
        repeat_tract_criteria: criteria.clone(),
    };
    let repeat_ground = match setting {
        RepeatGroundSetting::StrPaths => RepeatGround::StrPaths(sorting),
        // **The satellites are handed to the generic path as the segmentation is built**
        // ([`Segmentation::of_segments`]), so no route can record this setting without them.
        RepeatGroundSetting::CallSatellites => RepeatGround::CallSatellites(sorting),
        RepeatGroundSetting::GenericPathOnly => panic!(
            "a run that sorts no repeats is given its segments by segments_without_a_catalog, \
             not cut from a catalog"
        ),
    };
    Ok(Segmentation::of_segments(
        segments,
        analysed.clone(),
        repeat_ground,
    ))
}

/// Refuse a catalog that is not on disk, naming the reference it should have been built from.
///
/// **A file that is not there is not a file that will not read**: the second means a corrupt or
/// foreign catalog and sends the user to rebuild it, where this one means `repeat-catalog` has
/// not been run yet, or was run somewhere else.
fn refuse_a_catalog_that_is_not_there(catalog: &Path, reference: &Path) -> Result<(), GroundError> {
    match catalog.exists() {
        true => Ok(()),
        false => Err(GroundError::MissingCatalog {
            path: catalog.to_path_buf(),
            reference: reference.to_path_buf(),
        }),
    }
}

/// Render a catalog failure, naming the flag to move when the failure is that this run asked
/// for more than the file holds.
///
/// **The refusal is real and it is not policy** — the rows below the file's floors were never
/// written, so the request cannot be served (`run_ssr_observations.md` §2.3) — but on its own
/// it names two numbers and no way to change either. A person who typed `--min-copies 3,3,3,3,3,3`
/// should be told that flag, not left to infer which of six knobs produced *"period 1: catalog
/// holds tracts of 5 copies and up, reader asked for 3"*.
///
/// **The flank floor has no flag and so has no arm here**: a run pins it at the catalog's
/// own, so a catalog built at a wider flank than 15 bp is a file that has to be rebuilt, which
/// is what the general catalog error already says.
fn catalog_error_naming_the_flag(source: RepeatCatalogError, path: &Path) -> GroundError {
    // Exhaustive on the refusal on purpose: a new bounded axis must not silently inherit the
    // no-flag answer and leave a person hunting six knobs for the one they moved.
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

    /// A contig table of three contigs, as a psp header would hand it over.
    fn three_contigs() -> ContigList {
        ContigList {
            entries: [("chr1", 1_000), ("chr2", 500), ("chr3", 300)]
                .into_iter()
                .map(|(name, length)| crate::fasta::ContigEntry {
                    name: name.to_owned(),
                    length,
                    md5: None,
                })
                .collect(),
        }
    }

    /// The psps' ground: `chr1:101-500` and `chr2:1-500`, nothing on `chr3`.
    fn psps_ground(contigs: &ContigList) -> GenomeRegions {
        let directory = tempfile::tempdir().expect("a temporary directory");
        let bed = directory.path().join("walked.bed");
        std::fs::write(&bed, "chr1\t100\t500\nchr2\t0\t500\n").expect("the BED is written");
        let bounds = addressable_bounds(Path::new("first.psp"), contigs).expect("addressable");
        GenomeRegions::from_bed_path(&bed, &bounds).expect("the walked ground reads")
    }

    /// `--regions` given as `bed_text`, checked against [`psps_ground`].
    fn checked(bed_text: &str) -> Result<GenomeRegions, GroundError> {
        let contigs = three_contigs();
        let ground = psps_ground(&contigs);
        let directory = tempfile::tempdir().expect("a temporary directory");
        let bed = directory.path().join("asked.bed");
        std::fs::write(&bed, bed_text).expect("the BED is written");
        regions_inside_the_psps(&bed, Path::new("first.psp"), &contigs, &ground)
    }

    /// **Regions inside the psps' ground are accepted as the BED reader reads them**: lines out of
    /// order, overlapping and touching are merged, and a line running past its contig's end is cut
    /// there — so a BED written as overlapping pieces of the whole ground is the whole ground.
    #[test]
    fn regions_inside_the_psps_ground_are_accepted_merged() {
        let accepted = checked("chr2\t0\t9999\nchr1\t300\t500\nchr1\t100\t200\nchr1\t150\t300\n")
            .expect("every base lies in the psps' ground");
        let contigs = three_contigs();
        assert_eq!(accepted, psps_ground(&contigs));

        let inside = checked("chr1\t200\t210\n").expect("a piece inside");
        assert_eq!(inside.len(), 1);
    }

    /// **A BED reaching outside the psps' ground is refused with every piece outside named, in
    /// 1-based inclusive coordinates, beside what the psps hold on that contig** — the message a
    /// person reads, whole.
    #[test]
    fn regions_outside_the_psps_ground_are_refused_naming_each_piece() {
        let refused = checked("chr1\t50\t150\nchr1\t450\t600\nchr3\t0\t300\n")
            .expect_err("three pieces lie outside");

        let GroundError::RegionsOutsideThePsps {
            pieces, unlisted, ..
        } = &refused
        else {
            panic!("expected the regions-outside refusal, got {refused:?}");
        };
        assert_eq!(*unlisted, 0);
        assert_eq!(
            pieces
                .iter()
                .map(|piece| piece.piece.as_str())
                .collect::<Vec<_>>(),
            vec!["chr1:51-100", "chr1:501-600", "chr3:1-300"],
        );
        let rendered = refused.to_string();
        let expected = "\
--regions DIR/asked.bed asks for regions the psps do not cover, so nothing was called:
  chr1:51-100    on chr1 the psps cover chr1:101-500
  chr1:501-600   on chr1 the psps cover chr1:101-500
  chr3:1-300     the psps cover nothing on chr3
The psps cover 2 intervals, chr1:101-500 … chr2:1-500 (900 bases). Positions are \
1-based, both ends included.
Remove these regions from the BED, or write the psps again over regions that include them.";
        let bed_path = match &refused {
            GroundError::RegionsOutsideThePsps { bed, .. } => bed.parent().expect("a directory"),
            _ => unreachable!(),
        };
        assert_eq!(
            rendered.replace(&bed_path.display().to_string(), "DIR"),
            expected
        );
    }

    /// **A span bridging a gap in the psps' ground is refused for the gap**, and a contig with
    /// several walked spans is described by its first and last.
    #[test]
    fn a_region_bridging_a_gap_in_the_psps_ground_is_refused_for_the_gap() {
        let contigs = three_contigs();
        let directory = tempfile::tempdir().expect("a temporary directory");
        let walked = directory.path().join("walked.bed");
        std::fs::write(&walked, "chr1\t0\t100\nchr1\t200\t300\n").expect("written");
        let bounds = addressable_bounds(Path::new("first.psp"), &contigs).expect("addressable");
        let ground = GenomeRegions::from_bed_path(&walked, &bounds).expect("reads");
        let asked = directory.path().join("asked.bed");
        std::fs::write(&asked, "chr1\t50\t250\n").expect("written");

        let refused = regions_inside_the_psps(&asked, Path::new("first.psp"), &contigs, &ground)
            .expect_err("the gap is outside");

        let GroundError::RegionsOutsideThePsps { pieces, .. } = &refused else {
            panic!("expected the regions-outside refusal, got {refused:?}");
        };
        assert_eq!(
            pieces,
            &vec![RegionOutsideThePsps {
                piece: "chr1:101-200".to_owned(),
                what_the_psps_hold_there:
                    "on chr1 the psps cover 2 intervals, chr1:1-100 … chr1:201-300".to_owned(),
            }],
        );
    }

    /// **More than ten pieces outside are listed ten, then counted**, so a BED of thousands of
    /// lines does not fill a terminal.
    #[test]
    fn pieces_past_the_tenth_are_counted_not_listed() {
        let bed: String = (0..13)
            .map(|piece| format!("chr3\t{}\t{}\n", piece * 20, piece * 20 + 10))
            .collect();

        let refused = checked(&bed).expect_err("every piece is on a contig the psps lack");

        let GroundError::RegionsOutsideThePsps {
            pieces, unlisted, ..
        } = &refused
        else {
            panic!("expected the regions-outside refusal, got {refused:?}");
        };
        assert_eq!(pieces.len(), PIECES_LISTED);
        assert_eq!(*unlisted, 3);
        assert!(
            refused.to_string().contains("\n  … and 3 more\n"),
            "{refused}"
        );
    }

    /// **A contig the psps' table — and so the reference — lacks is the BED reader's refusal**, as
    /// for every other `--regions`, and so is a BED with no interval in it.
    #[test]
    fn a_bed_that_will_not_read_is_the_bed_readers_refusal() {
        let unknown = checked("chr1\t100\t200\nchrUn\t0\t10\n").expect_err("no such contig");
        let rendered = crate::error_render::format_error_chain(&unknown);
        assert!(
            matches!(
                unknown,
                GroundError::Bed {
                    source: BedError::UnknownContig { line_number: 2, .. },
                    ..
                }
            ),
            "{rendered}"
        );
        assert!(
            rendered.contains("contig `chrUn` is not in the reference contig list"),
            "{rendered}"
        );

        for empty in ["", "# only a comment\ntrack name=x\n"] {
            let refused = checked(empty).expect_err("an empty BED asks for nothing");
            assert!(
                matches!(
                    refused,
                    GroundError::Bed {
                        source: BedError::NoRegions,
                        ..
                    }
                ),
                "{refused:?}"
            );
        }
    }

    /// **A run asking for repeats the file does not hold is told which flag to move.**
    ///
    /// The catalog's own refusal names two numbers and no knob: *"period 1: catalog holds
    /// tracts of 5 copies and up, reader asked for 3"* leaves a person to work out which of
    /// six flags produced it. Each bounded axis maps to the flag that moves it.
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
