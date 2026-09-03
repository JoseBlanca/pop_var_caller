//! The `[segmentation]` section of the psp header: the ground the sample was analysed
//! over, and what the segmentation that shaped its observations was computed from.
//!
//! These are the fields a calling run's compatibility checks compare
//! (`doc/devel/ng/spec/run_streaming.md` §6.1, §6.2): the **analysed regions** are the
//! one cross-cohort check — two samples analysed over different ground are not
//! comparable — and the **repeat catalog's identity** plus the **repeat-tract routing
//! criteria** are the file-against-run check, because "no observation crosses a
//! segment's edge" holds only when the writing run's segmentation and the calling run's
//! are the same. The typed operand of both checks is
//! [`SegmentationInputs`], whose
//! [`first_difference`](SegmentationInputs::first_difference) names the field a refusal
//! reports — which is why this section records the inputs whole rather than a digest:
//! a digest can say two things disagree, never *what* disagrees.
//!
//! **Its own file, one seam of the header.** The wire types here are the TOML shape of
//! one header section, converted to and from the typed [`SegmentationInputs`]; the
//! framing, the version and every other section stay in
//! [`header`](crate::ng::psp::header), which is the only module that assembles a whole
//! header. Broken rules travel back as the header's own
//! [`BrokenRule`], so the writer refuses to produce what the reader would refuse to
//! believe, exactly as for every other field.

use serde::{Deserialize, Serialize};

use crate::ng::psp::header::{BrokenRule, ContigIdentity, MAX_TOML_INTEGER, digest_of, hex_of};
use crate::ng::reference_info::ContigInfo;
use crate::ng::region_typing::GenomeRegions;
use crate::ng::region_typing::segment_criteria::{MinCopies, SsrSegmentCriteria};
use crate::ng::repeat_catalog::{RepeatCatalogHeader, StrRepeatCriteria};
use crate::ng::run::segments::SegmentationInputs;
use crate::ng::tandem_repeat::{PeriodRange, ScanParams};
use crate::ng::types::{Bp, MAX_MOTIF_LEN};
use crate::regions::{ContigBounds, Region};

// ---------------------------------------------------------------------
// The wire shape — what the TOML body's [segmentation] section is
// ---------------------------------------------------------------------

/// The `[segmentation]` section, in TOML's own terms.
///
/// The same split as the header's other wire types: the typed
/// [`SegmentationInputs`] carries checked constructions — a [`PeriodRange`] that cannot
/// be empty, a [`GenomeRegions`] that cannot be out of order — and this carries what
/// TOML has, rebuilt through those checked constructors on the way in.
#[derive(Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub(crate) struct WireSegmentation {
    /// `[[segmentation.analysed-region]]` — the ground, in genomic order, each span
    /// naming its contig by the header's own contig list.
    #[serde(default)]
    analysed_region: Vec<WireAnalysedRegion>,
    /// `[segmentation.repeat-tract-criteria]` — the criteria this run asked the
    /// catalog with when routing the ground into segments, which decide where a
    /// segment ends (the same value [`SegmentationInputs::repeat_tract_criteria`]
    /// holds in memory).
    repeat_tract_criteria: WireRepeatCriteria,
    /// `[segmentation.catalog]` — the repeat catalog's own header, whole, so a refusal
    /// can name the exact field of it that differs.
    catalog: WireCatalog,
}

/// One analysed span: 1-based inclusive, the same convention as every region in the
/// pipeline (`src/regions.rs`).
#[derive(Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
struct WireAnalysedRegion {
    /// The contig's name — resolved against the header's own contig list on the way in,
    /// so a span cannot quietly index a contig the file does not declare.
    contig: String,
    /// 1-based inclusive lower bound.
    start: u64,
    /// 1-based inclusive upper bound.
    end: u64,
}

/// A [`StrRepeatCriteria`] in TOML's terms — used twice: once for the criteria the run
/// routed with, once for the criteria the catalog records it was built under.
#[derive(Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
struct WireRepeatCriteria {
    /// The smallest repeat period classified, inclusive (`>= 1`).
    period_min: u8,
    /// The largest repeat period classified, inclusive (`>= period-min`, at most
    /// [`MAX_MOTIF_LEN`]).
    period_max: u8,
    /// The fewest motif copies a tract needs, per period: entry `p - 1` is the floor
    /// for period `p`. Exactly [`MAX_MOTIF_LEN`] entries.
    min_copies_by_period: Vec<u32>,
    /// The copy floor for a period past the table's end.
    min_copies_for_wider_periods: u32,
    /// The purity floor, a fraction in `[0, 1]` — enforced by the rules, not only
    /// documented.
    min_purity: f64,
    /// The detector-score floor.
    min_score: i32,
    /// The bundle-clustering radius, in bases (`>= 1`; a radius of zero bundles
    /// nothing and the classifier asserts against it).
    bundle_threshold_bp: u64,
    /// Sequence required on each side of a tract, in bases.
    min_flank_bp: u64,
    /// A tract longer than this is a satellite, not a locus, in bases.
    max_str_len_bp: u64,
}

/// A [`RepeatCatalogHeader`] in TOML's terms. Every field of it, because the operand of
/// the file-against-run check is whole-header equality and the refusal names the first
/// field that differs.
#[derive(Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
struct WireCatalog {
    /// The whole-reference digest the catalog was built on, 32 lowercase hex.
    reference_md5: String,
    /// The crate version that wrote the catalog.
    tool_version: String,
    /// `[segmentation.catalog.scan]` — the scoring weights every stored tract came out
    /// of.
    scan: WireScan,
    /// `[segmentation.catalog.built-under]` — the criteria the catalog's builder was
    /// given.
    built_under: WireRepeatCriteria,
    /// `[[segmentation.catalog.contig]]` — the catalog's own contig table, geometry and
    /// all, with each contig's longest stored tract beside it. **The catalog's record,
    /// not this file's**: the psp's own coordinate space is the header's top-level
    /// contig list, and the two are compared by the run, not conflated here.
    #[serde(default)]
    contig: Vec<WireCatalogContig>,
}

/// The lag-`p` scoring weights ([`ScanParams`]).
#[derive(Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
struct WireScan {
    match_reward: i32,
    mismatch_penalty: i32,
    min_copies: u32,
}

/// One contig as the catalog records it ([`ContigInfo`]), plus its longest stored
/// tract — beside the contig rather than in a parallel list, so the two cannot go
/// different lengths in the file.
#[derive(Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
struct WireCatalogContig {
    name: String,
    length: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    md5: Option<String>,
    /// `.fai` geometry of the reference the catalog was built from — part of the
    /// catalog header's identity, compared whole.
    offset: u64,
    line_bases: u64,
    line_width: u64,
    /// The longest tract stored in this contig, in bases; 0 where it holds none.
    longest_tract_bp: u64,
}

// ---------------------------------------------------------------------
// Typed → wire
// ---------------------------------------------------------------------

impl WireSegmentation {
    /// The wire shape of `inputs`, spelling each analysed span's contig from the
    /// header's own contig list.
    ///
    /// The writer's rules ([`check_segmentation`]) refuse a span outside that list
    /// before an honest encode reaches here. The fallback spelling below exists for the
    /// tests that deliberately serialise a rule-breaking header to meet it from the
    /// reading side: no `encode`-produced file can carry it, and a hand-built file that
    /// does is refused on the way back in rather than quietly re-anchored.
    pub(crate) fn from_inputs(inputs: &SegmentationInputs, contigs: &[ContigIdentity]) -> Self {
        // Exhaustive, so a field added to the operand type in a later step fails to
        // compile *here*, on the encode side — not only in the decode literal, where a
        // default could be filled in without the writer ever recording the field.
        let SegmentationInputs {
            catalog,
            repeat_tract_criteria,
            analysed_regions,
        } = inputs;
        let analysed_region = analysed_regions
            .spans()
            .iter()
            .map(|span| {
                let index = span.chrom_id as usize;
                let contig = contigs
                    .get(index)
                    .map(|contig| contig.name.clone())
                    .unwrap_or_else(|| format!("undeclared-contig-{index}"));
                WireAnalysedRegion {
                    contig,
                    start: u64::from(span.start),
                    end: u64::from(span.end),
                }
            })
            .collect();
        WireSegmentation {
            analysed_region,
            repeat_tract_criteria: WireRepeatCriteria::from_criteria(repeat_tract_criteria),
            catalog: WireCatalog::from_catalog(catalog),
        }
    }
}

impl WireRepeatCriteria {
    fn from_criteria(criteria: &StrRepeatCriteria) -> Self {
        // Exhaustive for the same reason as `from_inputs`: a new criterion must reach
        // the wire or fail to compile.
        let StrRepeatCriteria {
            classification,
            min_flank_bp,
            max_str_len_bp,
        } = criteria;
        let SsrSegmentCriteria {
            periods,
            min_copies,
            min_purity,
            min_score,
            bundle_threshold,
        } = classification;
        WireRepeatCriteria {
            period_min: periods.min(),
            period_max: periods.max(),
            min_copies_by_period: (1..=MAX_MOTIF_LEN as u8)
                .map(|period| min_copies.for_period(period))
                .collect(),
            min_copies_for_wider_periods: min_copies.for_wider_periods(),
            min_purity: wire_float_of(*min_purity),
            min_score: *min_score,
            bundle_threshold_bp: *bundle_threshold,
            min_flank_bp: min_flank_bp.get(),
            max_str_len_bp: max_str_len_bp.get(),
        }
    }
}

impl WireCatalog {
    fn from_catalog(catalog: &RepeatCatalogHeader) -> Self {
        // Exhaustive for the same reason as `from_inputs`.
        let RepeatCatalogHeader {
            contigs,
            reference_md5,
            built_under,
            scan,
            tool_version,
            longest_tract_bp,
        } = catalog;
        WireCatalog {
            reference_md5: hex_of(*reference_md5),
            tool_version: tool_version.clone(),
            scan: WireScan {
                match_reward: scan.match_reward,
                mismatch_penalty: scan.mismatch_penalty,
                min_copies: scan.min_copies,
            },
            built_under: WireRepeatCriteria::from_criteria(built_under),
            contig: contigs
                .iter()
                .enumerate()
                .map(|(index, contig)| WireCatalogContig {
                    name: contig.name.clone(),
                    length: contig.length,
                    md5: contig.md5.map(hex_of),
                    offset: contig.offset,
                    line_bases: contig.line_bases,
                    line_width: contig.line_width,
                    // The rules refuse unequal lists before an honest encode reaches
                    // here; the fallback keeps this callable on a deliberately broken
                    // header, whose file the rules then refuse from the reading side.
                    longest_tract_bp: longest_tract_bp.get(index).copied().unwrap_or(0),
                })
                .collect(),
        }
    }
}

/// An `f32` widened for TOML through its own shortest decimal, so the header shows
/// `0.93` rather than `0.9300000071525574`.
///
/// Exact both ways: `Display` on an `f32` prints the shortest decimal that reads back
/// to the same `f32`, and that decimal's nearest `f64` narrows back to it.
///
/// The `expect` cannot fire on any input: `f64`'s parser accepts every string `f32`'s
/// `Display` produces — `NaN` and `inf` included — so even the test-only path that
/// serialises a rule-breaking header (the rules refuse non-finite values on every
/// honest encode) parses back here.
fn wire_float_of(value: f32) -> f64 {
    format!("{value}")
        .parse()
        .expect("a float's own Display re-parses")
}

// ---------------------------------------------------------------------
// Wire → typed
// ---------------------------------------------------------------------

impl WireSegmentation {
    /// The typed inputs, rebuilt through the same checked constructors the run itself
    /// uses — a period range that cannot be empty, a region set that cannot be out of
    /// order — with `contigs` (the header's own list) anchoring each span's contig
    /// name back to a coordinate.
    pub(crate) fn into_inputs(
        self,
        contigs: &[ContigIdentity],
    ) -> Result<SegmentationInputs, BrokenRule> {
        let spans = self
            .analysed_region
            .iter()
            .map(|region| region.to_span(contigs))
            .collect::<Result<Vec<_>, BrokenRule>>()?;
        let bounds = contig_bounds_of(contigs);
        let analysed_regions =
            GenomeRegions::from_normalized_spans(spans, &bounds).map_err(|broken| {
                BrokenRule::new("segmentation.analysed-region", broken.to_string())
            })?;

        Ok(SegmentationInputs {
            catalog: self.catalog.into_catalog()?,
            repeat_tract_criteria: self
                .repeat_tract_criteria
                .into_criteria("segmentation.repeat-tract-criteria")?,
            analysed_regions,
        })
    }
}

/// The header's contig list in the shape the span constructors take.
///
/// The pipeline's spans are 32-bit; a longer contig cannot hold a recordable span past
/// `u32::MAX` anyway ([`WireAnalysedRegion::to_span`] refused it), so **clamping** the
/// bound loses no check — where a *truncating* cast would turn a >4 Gbp contig's bound
/// into a small number and refuse every span the writer accepted.
pub(crate) fn contig_bounds_of(contigs: &[ContigIdentity]) -> Vec<ContigBounds<'_>> {
    contigs
        .iter()
        .map(|contig| ContigBounds {
            name: &contig.name,
            length: contig.length.min(u64::from(u32::MAX)) as u32,
        })
        .collect()
}

impl WireAnalysedRegion {
    /// This span in the pipeline's own terms, its contig name resolved to its position
    /// in the header's contig list.
    fn to_span(&self, contigs: &[ContigIdentity]) -> Result<Region, BrokenRule> {
        let field = "segmentation.analysed-region";
        let Some(chrom_id) = contigs.iter().position(|contig| contig.name == self.contig) else {
            return Err(BrokenRule::new(
                field,
                format!(
                    "{:?} is not in this file's contig list; an analysed span must name \
                     one of the file's own contigs",
                    self.contig
                ),
            ));
        };
        let coordinate = |name: &str, value: u64| {
            u32::try_from(value).map_err(|_| {
                BrokenRule::new(
                    field,
                    format!(
                        "{name} {value} on {:?} is wider than the 32-bit coordinates \
                         analysed spans are stored in",
                        self.contig
                    ),
                )
            })
        };
        Ok(Region {
            // A file of more than `u32::MAX` contigs cannot exist: its header would
            // exceed the body ceiling long before.
            chrom_id: chrom_id as u32,
            start: coordinate("start", self.start)?,
            end: coordinate("end", self.end)?,
        })
    }
}

impl WireRepeatCriteria {
    fn into_criteria(self, section: &str) -> Result<StrRepeatCriteria, BrokenRule> {
        let periods = PeriodRange::new(self.period_min, self.period_max).map_err(|refused| {
            BrokenRule::new(format!("{section}.period-min"), refused.to_string())
        })?;
        let by_period: [u32; MAX_MOTIF_LEN] =
            self.min_copies_by_period
                .try_into()
                .map_err(|found: Vec<u32>| {
                    BrokenRule::new(
                        format!("{section}.min-copies-by-period"),
                        format!(
                            "has {} entries; the table is {MAX_MOTIF_LEN} wide",
                            found.len()
                        ),
                    )
                })?;
        Ok(StrRepeatCriteria {
            classification: SsrSegmentCriteria {
                periods,
                min_copies: MinCopies::new(by_period, self.min_copies_for_wider_periods),
                // Exact for every value this build wrote (`wire_float_of`); a hand-edited
                // value narrows to the nearest `f32`, and an out-of-range one is refused
                // by the rules right after this conversion.
                min_purity: self.min_purity as f32,
                min_score: self.min_score,
                bundle_threshold: self.bundle_threshold_bp,
            },
            min_flank_bp: Bp(self.min_flank_bp),
            max_str_len_bp: Bp(self.max_str_len_bp),
        })
    }
}

impl WireCatalog {
    fn into_catalog(self) -> Result<RepeatCatalogHeader, BrokenRule> {
        let mut contigs = Vec::with_capacity(self.contig.len());
        let mut longest_tract_bp = Vec::with_capacity(self.contig.len());
        for contig in self.contig {
            // The digest is checked while the name is still at hand, so a refusal can
            // say *which* catalog row broke — at 30,000 scaffolds "a digest is wrong
            // somewhere" is not actionable.
            let md5 = contig
                .md5
                .map(|spelled| {
                    digest_of(
                        &format!("segmentation.catalog.contig.{}.md5", contig.name),
                        &spelled,
                    )
                })
                .transpose()?;
            contigs.push(ContigInfo {
                md5,
                name: contig.name,
                length: contig.length,
                offset: contig.offset,
                line_bases: contig.line_bases,
                line_width: contig.line_width,
            });
            longest_tract_bp.push(contig.longest_tract_bp);
        }
        Ok(RepeatCatalogHeader {
            contigs,
            longest_tract_bp,
            reference_md5: digest_of("segmentation.catalog.reference-md5", &self.reference_md5)?,
            built_under: self
                .built_under
                .into_criteria("segmentation.catalog.built-under")?,
            scan: ScanParams {
                match_reward: self.scan.match_reward,
                mismatch_penalty: self.scan.mismatch_penalty,
                min_copies: self.scan.min_copies,
            },
            tool_version: self.tool_version,
        })
    }
}

// ---------------------------------------------------------------------
// The rules — checked on both sides, like every other header rule
// ---------------------------------------------------------------------

/// Every rule the segmentation section must satisfy, called from the header's own
/// `check_rules` so writer and reader refuse the same things.
pub(crate) fn check_segmentation(
    inputs: &SegmentationInputs,
    contigs: &[ContigIdentity],
) -> Result<(), BrokenRule> {
    if inputs.analysed_regions.is_empty() {
        return Err(BrokenRule::new(
            "segmentation.analysed-region",
            "is empty; a psp records the ground its sample was analysed over",
        ));
    }
    for span in inputs.analysed_regions.spans() {
        let Some(contig) = contigs.get(span.chrom_id as usize) else {
            return Err(BrokenRule::new(
                "segmentation.analysed-region",
                format!(
                    "a span indexes contig {} and this file declares {}; the analysed \
                     regions must be built against the file's own contig list",
                    span.chrom_id,
                    contigs.len()
                ),
            ));
        };
        if u64::from(span.end) > contig.length {
            return Err(BrokenRule::new(
                "segmentation.analysed-region",
                format!(
                    "a span ends at {} on {}, which is {} bases long",
                    span.end, contig.name, contig.length
                ),
            ));
        }
    }

    check_criteria(
        "segmentation.repeat-tract-criteria",
        &inputs.repeat_tract_criteria,
    )?;
    check_catalog(&inputs.catalog)
}

/// The criteria's own invariants, as the classifier release-asserts them
/// (`segment_criteria.rs::classify`) — enforced here so a file carrying a value the
/// pipeline could never produce is refused naming its field, instead of surviving as a
/// comparison operand until some later consumer panics on it.
fn check_criteria(section: &str, criteria: &StrRepeatCriteria) -> Result<(), BrokenRule> {
    let purity = criteria.classification.min_purity;
    if !(0.0..=1.0).contains(&purity) {
        return Err(BrokenRule::new(
            format!("{section}.min-purity"),
            format!("is {purity}; a purity floor is a fraction between 0 and 1"),
        ));
    }
    if criteria.classification.bundle_threshold < 1 {
        return Err(BrokenRule::new(
            format!("{section}.bundle-threshold-bp"),
            "is 0; a bundle radius of zero bundles nothing and the classifier refuses it",
        ));
    }
    if usize::from(criteria.classification.periods.max()) > MAX_MOTIF_LEN {
        return Err(BrokenRule::new(
            format!("{section}.period-max"),
            format!(
                "is {}; a repeat motif is at most {MAX_MOTIF_LEN} bases",
                criteria.classification.periods.max()
            ),
        ));
    }
    check_toml_integer(
        &format!("{section}.bundle-threshold-bp"),
        criteria.classification.bundle_threshold,
    )?;
    check_toml_integer(
        &format!("{section}.min-flank-bp"),
        criteria.min_flank_bp.get(),
    )?;
    check_toml_integer(
        &format!("{section}.max-str-len-bp"),
        criteria.max_str_len_bp.get(),
    )
}

fn check_catalog(catalog: &RepeatCatalogHeader) -> Result<(), BrokenRule> {
    // A newline in one of the catalog's free strings would land in the body as further
    // lines — the forged-key hole the header's field-name and contig-name rules
    // document having closed once already. The catalog's rows are writer-controlled
    // strings the same way, so they get the same rule.
    check_plain_single_line("segmentation.catalog.tool-version", &catalog.tool_version)?;
    if catalog.contigs.len() != catalog.longest_tract_bp.len() {
        return Err(BrokenRule::new(
            "segmentation.catalog.contig",
            format!(
                "the catalog names {} contigs and {} longest-tract lengths; the two \
                 lists travel row by row and must match",
                catalog.contigs.len(),
                catalog.longest_tract_bp.len()
            ),
        ));
    }
    for (contig, longest) in catalog.contigs.iter().zip(&catalog.longest_tract_bp) {
        let field = |name: &str| format!("segmentation.catalog.contig.{name}");
        check_plain_single_line(&field("name"), &contig.name)?;
        check_toml_integer(&field("length"), contig.length)?;
        check_toml_integer(&field("offset"), contig.offset)?;
        check_toml_integer(&field("line-bases"), contig.line_bases)?;
        check_toml_integer(&field("line-width"), contig.line_width)?;
        check_toml_integer(&field("longest-tract-bp"), *longest)?;
    }
    check_criteria("segmentation.catalog.built-under", &catalog.built_under)
}

/// A writer-controlled string that appears in the header's text may not carry
/// whitespace or a control character: TOML writes such a value as a multi-line string,
/// so its own bytes land in the body as further lines and can make `head` show a key no
/// field declared. The same rule, for the same reason, as the header's contig-name and
/// manifest-field-name checks.
fn check_plain_single_line(field: &str, value: &str) -> Result<(), BrokenRule> {
    if value
        .chars()
        .any(|character| character.is_whitespace() || character.is_control())
    {
        return Err(BrokenRule::new(
            field.to_string(),
            format!(
                "{value:?} holds whitespace or a control character; such a value is \
                 written as a multi-line string and can make the header's text show a \
                 key no field declares"
            ),
        ));
    }
    Ok(())
}

/// TOML's integer is signed 64-bit, so a `u64` above `i64::MAX` would serialize into a
/// file the writer's own reader refuses — the same rule the header applies to its
/// contig lengths.
fn check_toml_integer(field: &str, value: u64) -> Result<(), BrokenRule> {
    if value > MAX_TOML_INTEGER {
        return Err(BrokenRule::new(
            field.to_string(),
            format!(
                "is {value}; a TOML integer is signed, so a header cannot carry more \
                 than {MAX_TOML_INTEGER}"
            ),
        ));
    }
    Ok(())
}

// ---------------------------------------------------------------------
// Fixtures — shared with the header's and the writer's tests
// ---------------------------------------------------------------------

/// A [`SegmentationInputs`] for tests, anchored to a header's contig list.
///
/// **Every value differs from its type's default wherever the type has one**, so a
/// round trip that quietly substituted a default fails rather than passing by
/// coincidence: the analysed regions are proper sub-spans on any contig four bases or
/// longer, the routing purity floor is 0.93, the catalog's build criteria carry a
/// non-default flank floor and score floor, its scan weights differ in all three
/// fields, and its builder version is its own.
#[cfg(test)]
pub(crate) fn segmentation_inputs_for_tests(contigs: &[ContigIdentity]) -> SegmentationInputs {
    let bounds = contig_bounds_of(contigs);
    let spans = contigs
        .iter()
        .zip(&bounds)
        .enumerate()
        .map(|(index, (_, bound))| {
            // A quarter in from each end, so the span is a proper sub-span of any
            // contig long enough to have one and stays valid down to one base.
            let quarter = bound.length / 4;
            let start = quarter.max(1);
            Region {
                chrom_id: index as u32,
                start,
                end: (bound.length - quarter).max(start),
            }
        })
        .collect();
    let analysed_regions = GenomeRegions::from_normalized_spans(spans, &bounds)
        .expect("the fixture's spans are normalized");

    let mut repeat_tract_criteria = StrRepeatCriteria::default();
    repeat_tract_criteria.classification.min_purity = 0.93;

    let mut built_under = StrRepeatCriteria {
        min_flank_bp: Bp(11),
        ..StrRepeatCriteria::default()
    };
    built_under.classification.min_score = 9;

    let catalog_contigs: Vec<ContigInfo> = contigs
        .iter()
        .enumerate()
        .map(|(index, contig)| ContigInfo {
            name: contig.name.clone(),
            length: contig.length,
            offset: 60 + index as u64 * 7,
            line_bases: 60,
            line_width: 61,
            md5: contig.md5,
        })
        .collect();
    let longest_tract_bp = (0..catalog_contigs.len()).map(|i| 88 + i as u64).collect();
    SegmentationInputs {
        catalog: RepeatCatalogHeader {
            contigs: catalog_contigs,
            longest_tract_bp,
            reference_md5: [0x0a; 16],
            built_under,
            scan: ScanParams {
                match_reward: 3,
                mismatch_penalty: 5,
                min_copies: 4,
            },
            tool_version: "trf-port-0.9.1".to_string(),
        },
        repeat_tract_criteria,
        analysed_regions,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::tandem_repeat::{
        DEFAULT_MATCH_REWARD, DEFAULT_MIN_COPIES, DEFAULT_MISMATCH_PENALTY,
    };

    fn contigs() -> Vec<ContigIdentity> {
        vec![
            ContigIdentity {
                name: "SL4.0ch01".to_string(),
                length: 90_863_682,
                md5: Some([0x1b; 16]),
            },
            ContigIdentity {
                name: "SL4.0ch02".to_string(),
                length: 53_473_368,
                md5: Some([0x2c; 16]),
            },
        ]
    }

    /// The section's own round trip, field for field — the header's tests then prove
    /// the same thing through the whole encode/decode path.
    #[test]
    fn the_section_round_trips_field_for_field() {
        let contigs = contigs();
        let written = segmentation_inputs_for_tests(&contigs);
        let read = WireSegmentation::from_inputs(&written, &contigs)
            .into_inputs(&contigs)
            .expect("this build's own section reads back");
        assert_eq!(read, written);
    }

    /// The fixture keeps the promise its doc makes: the three operands of the future
    /// compatibility check — routing criteria, catalog build criteria, scan weights —
    /// all differ from their defaults, so a decode that substituted a default for any
    /// of them fails the round trip above instead of passing by coincidence.
    #[test]
    fn the_fixture_differs_from_every_default_the_operands_have() {
        let inputs = segmentation_inputs_for_tests(&contigs());
        assert_ne!(inputs.repeat_tract_criteria, StrRepeatCriteria::default());
        assert_ne!(inputs.catalog.built_under, StrRepeatCriteria::default());
        assert_ne!(inputs.catalog.scan, ScanParams::default());
        // All three scan weights individually, so a partial substitution cannot hide
        // behind the struct-level inequality.
        assert_ne!(inputs.catalog.scan.match_reward, DEFAULT_MATCH_REWARD);
        assert_ne!(
            inputs.catalog.scan.mismatch_penalty,
            DEFAULT_MISMATCH_PENALTY
        );
        assert_ne!(inputs.catalog.scan.min_copies, DEFAULT_MIN_COPIES);
    }

    /// The purity floor survives the trip through TOML's `f64` bit-for-bit, and its
    /// wire spelling is the short decimal a person wrote, not the widened one.
    #[test]
    fn a_purity_floor_stays_exact_and_short_through_the_wire_float() {
        let widened = wire_float_of(0.93_f32);
        assert_eq!(format!("{widened}"), "0.93");
        assert_eq!(widened as f32, 0.93_f32);
    }

    /// The whole-contig span — every whole-genome run's shape — passes the rules and
    /// round-trips. The shared fixture deliberately uses proper sub-spans, so without
    /// this test a `>` → `>=` regression in the span-end rule would refuse every
    /// default run while the suite stayed green.
    #[test]
    fn a_whole_contig_analysed_span_round_trips() {
        let contigs = contigs();
        let mut inputs = segmentation_inputs_for_tests(&contigs);
        inputs.analysed_regions = GenomeRegions::whole_contigs(&contig_bounds_of(&contigs));
        check_segmentation(&inputs, &contigs).expect("a whole-genome ground is legal");
        let read = WireSegmentation::from_inputs(&inputs, &contigs)
            .into_inputs(&contigs)
            .expect("and reads back");
        assert_eq!(read, inputs);
    }

    /// A contig longer than the 32-bit span space still anchors its (32-bit) spans:
    /// the bound is clamped, not truncated — at 2^32 a truncating cast gives a
    /// 0-length bound and refuses every span the writer accepted.
    #[test]
    fn a_span_on_a_contig_longer_than_u32_max_round_trips() {
        let contigs = vec![ContigIdentity {
            name: "huge".to_string(),
            length: 1_u64 << 32,
            md5: None,
        }];
        let mut inputs = segmentation_inputs_for_tests(&contigs);
        let bounds = [ContigBounds {
            name: "huge",
            length: u32::MAX,
        }];
        inputs.analysed_regions = GenomeRegions::from_normalized_spans(
            vec![Region {
                chrom_id: 0,
                start: 100,
                end: 1_000,
            }],
            &bounds,
        )
        .expect("a normalized span");
        check_segmentation(&inputs, &contigs).expect("the writer accepts it");
        let read = WireSegmentation::from_inputs(&inputs, &contigs)
            .into_inputs(&contigs)
            .expect("and the reader accepts the same file");
        assert_eq!(read, inputs);
    }

    /// A value at exactly the TOML integer ceiling is legal and one past it is not —
    /// the boundary itself, not only a value far beyond it, because a tightened bound
    /// refusing a legal file would otherwise pass the suite.
    #[test]
    fn a_criteria_value_at_the_toml_ceiling_is_accepted_and_one_more_refused() {
        let contigs = contigs();
        let mut inputs = segmentation_inputs_for_tests(&contigs);
        inputs.repeat_tract_criteria.min_flank_bp = Bp(MAX_TOML_INTEGER);
        check_segmentation(&inputs, &contigs).expect("the widest TOML integer is legal");
        let read = WireSegmentation::from_inputs(&inputs, &contigs)
            .into_inputs(&contigs)
            .expect("and round-trips");
        assert_eq!(read, inputs);

        inputs.repeat_tract_criteria.min_flank_bp = Bp(MAX_TOML_INTEGER + 1);
        let broken = check_segmentation(&inputs, &contigs).expect_err("one more is refused");
        assert_eq!(
            broken.field,
            "segmentation.repeat-tract-criteria.min-flank-bp"
        );
    }

    /// A span naming a contig the file does not declare is refused naming the span's
    /// own contig, not swallowed into a wrong coordinate.
    #[test]
    fn a_span_on_an_undeclared_contig_is_refused() {
        let contigs = contigs();
        let mut wire =
            WireSegmentation::from_inputs(&segmentation_inputs_for_tests(&contigs), &contigs);
        wire.analysed_region[0].contig = "SL4.0ch03".to_string();
        let refused = wire
            .into_inputs(&contigs)
            .expect_err("an undeclared contig is refused");
        assert!(refused.reason.contains("SL4.0ch03"), "{}", refused.reason);
    }

    /// A copy-floor table of the wrong width is refused; believing it would silently
    /// shift every period's floor.
    #[test]
    fn a_copy_floor_table_of_the_wrong_width_is_refused() {
        let contigs = contigs();
        let mut wire =
            WireSegmentation::from_inputs(&segmentation_inputs_for_tests(&contigs), &contigs);
        wire.repeat_tract_criteria.min_copies_by_period.pop();
        let refused = wire
            .into_inputs(&contigs)
            .expect_err("a short table is refused");
        assert!(
            refused.field.ends_with("min-copies-by-period"),
            "{}: {}",
            refused.field,
            refused.reason
        );
    }

    /// Both broken period ranges are refused through the same checked constructor the
    /// run itself builds with, and each refusal names the period field.
    #[test]
    fn a_broken_period_range_is_refused_naming_period_min() {
        let contigs = contigs();

        let mut inverted =
            WireSegmentation::from_inputs(&segmentation_inputs_for_tests(&contigs), &contigs);
        inverted.repeat_tract_criteria.period_min = 5;
        inverted.repeat_tract_criteria.period_max = 3;
        let refused = inverted
            .into_inputs(&contigs)
            .expect_err("an inverted range is refused");
        assert!(
            refused.field.ends_with("period-min"),
            "{}: {}",
            refused.field,
            refused.reason
        );

        let mut zero =
            WireSegmentation::from_inputs(&segmentation_inputs_for_tests(&contigs), &contigs);
        zero.repeat_tract_criteria.period_min = 0;
        let refused = zero
            .into_inputs(&contigs)
            .expect_err("a zero period floor is refused");
        assert!(
            refused.field.ends_with("period-min"),
            "{}: {}",
            refused.field,
            refused.reason
        );
    }

    /// The rules refuse what the constructors cannot represent wrongly, naming the
    /// field each time.
    #[test]
    fn the_rules_refuse_each_broken_shape_naming_its_field() {
        let contigs = contigs();
        let refused = |mangle: &dyn Fn(&mut SegmentationInputs)| {
            let mut inputs = segmentation_inputs_for_tests(&contigs);
            mangle(&mut inputs);
            check_segmentation(&inputs, &contigs).expect_err("the broken shape is refused")
        };

        // The classifier's own invariants, met from a file: a purity outside [0, 1]
        // (NaN included), a zero bundle radius, a period past the widest motif.
        let broken = refused(&|inputs| {
            inputs.repeat_tract_criteria.classification.min_purity = f32::NAN;
        });
        assert_eq!(
            broken.field,
            "segmentation.repeat-tract-criteria.min-purity"
        );

        let broken = refused(&|inputs| {
            inputs.repeat_tract_criteria.classification.min_purity = 7.0;
        });
        assert_eq!(
            broken.field,
            "segmentation.repeat-tract-criteria.min-purity"
        );

        let broken = refused(&|inputs| {
            inputs.repeat_tract_criteria.classification.bundle_threshold = 0;
        });
        assert_eq!(
            broken.field,
            "segmentation.repeat-tract-criteria.bundle-threshold-bp"
        );

        let broken = refused(&|inputs| {
            inputs.catalog.built_under.classification.periods =
                PeriodRange::new(1, MAX_MOTIF_LEN as u8 + 3).expect("a range the scanner allows");
        });
        assert_eq!(broken.field, "segmentation.catalog.built-under.period-max");

        let broken = refused(&|inputs| {
            inputs.catalog.longest_tract_bp.pop();
        });
        assert_eq!(broken.field, "segmentation.catalog.contig");

        let broken = refused(&|inputs| {
            inputs.catalog.built_under.max_str_len_bp = Bp(u64::MAX);
        });
        assert_eq!(
            broken.field,
            "segmentation.catalog.built-under.max-str-len-bp"
        );

        // The forged-line hole: a newline in a catalog contig name or the tool version
        // would land in the body as further lines.
        let broken = refused(&|inputs| {
            inputs.catalog.contigs[0].name = "evil\nfabricated-key = \"x\"".to_string();
        });
        assert_eq!(broken.field, "segmentation.catalog.contig.name");

        let broken = refused(&|inputs| {
            inputs.catalog.tool_version = "trf\nport".to_string();
        });
        assert_eq!(broken.field, "segmentation.catalog.tool-version");

        // Regions checked against a *different* contig list than they were built for:
        // the second contig's span now reaches past the shortened list's end.
        let with_a_shortened_second_contig = vec![contigs[0].clone(), {
            let mut shortened = contigs[1].clone();
            shortened.length = 1_000;
            shortened
        }];
        let inputs = segmentation_inputs_for_tests(&contigs);
        let broken = check_segmentation(&inputs, &with_a_shortened_second_contig)
            .expect_err("a span past its contig is refused");
        assert_eq!(broken.field, "segmentation.analysed-region");

        let broken = check_segmentation(&inputs, &contigs[..1])
            .expect_err("a span off the contig list is refused");
        assert_eq!(broken.field, "segmentation.analysed-region");
    }
}
