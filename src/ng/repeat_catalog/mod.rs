//! The reference's tandem-repeat catalog — every repeat in a genome, found once and
//! written beside the FASTA, so that a reader derives the genome's segmentation and its
//! STR loci under any copy floor it chooses without opening the reference again.
//!
//! It exists because the parameter pre-pass draws a **random sample of STR loci from each
//! (period, repeat count) stratum**, and "the `cap` lowest hashes in a stratum" needs the
//! stratum enumerated first — a genome-wide list nothing in ng could produce. Spec:
//! `doc/devel/ng/spec/repeat_catalog.md`; types and contracts:
//! `doc/devel/ng/arch/repeat_catalog.md`.
//!
//! **It is a catalog, and that is all it is.** Which of these tracts is an STR for calling
//! purposes is the reader's decision, made with [`StrRepeatCriteria`]; the file applies no
//! purity floor, no satellite cap and no bundling. The two thresholds it does apply — a
//! per-period copy floor of `[5, 5, 4, 4, 4, 3]` and 15 bp of sequence beside each tract —
//! sit below every floor a caller routes on, so every routing question stays answerable by
//! filtering the file (spec §1, §4.1).

pub mod builder;
pub mod criteria;
pub mod parquet_file;
pub mod reader;
pub mod row;
pub mod segments;
pub mod strata;

use std::path::PathBuf;

use crate::ng::reference_info::ContigInfo;
use crate::ng::tandem_repeat::ScanParams;
use crate::ng::types::{ContigId, Motif, Position};

pub use builder::{BuildTally, RepeatCatalogBuilder};
pub use criteria::{CriteriaRefusal, StrRepeatCriteria};
pub use parquet_file::RepeatCatalogWriter;
pub use reader::RepeatCatalog;
pub use row::{RowRejection, row_for_interval};
pub use strata::{StratumCounts, StratumSample};

/// How many rows a build wrote, per period — what `repeat-catalog` prints when it finishes
/// (spec §2.6), and the measurement open question 1 asks for.
///
/// Periods beyond [`MAX_MOTIF_LEN`](crate::ng::types::MAX_MOTIF_LEN) are counted together:
/// a scan cannot emit one, and a tally with a slot nothing can reach explains less than a
/// tally that says so.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct RowsByPeriod {
    by_period: [u64; crate::ng::types::MAX_MOTIF_LEN],
    beyond_table: u64,
}

impl RowsByPeriod {
    /// Count one row of `period`.
    #[inline]
    pub fn count(&mut self, period: u8) {
        match period
            .checked_sub(1)
            .and_then(|i| self.by_period.get_mut(usize::from(i)))
        {
            Some(slot) => *slot += 1,
            None => self.beyond_table += 1,
        }
    }

    /// Rows written at `period`.
    #[inline]
    pub fn for_period(&self, period: u8) -> u64 {
        period
            .checked_sub(1)
            .and_then(|i| self.by_period.get(usize::from(i)))
            .copied()
            .unwrap_or(self.beyond_table)
    }

    /// Rows written, all periods.
    #[inline]
    pub fn total(&self) -> u64 {
        self.by_period.iter().sum::<u64>() + self.beyond_table
    }
}

/// A tandem-repeat tract's span on one contig: **1-based inclusive**, ng's convention
/// (`ng_step_interfaces.md` §1).
///
/// A named pair rather than two bare `u64`s because a row carries **two** of them — the
/// span the scanner reported and the same tract cut back to whole motif copies — and the
/// two are adjacent, same-typed and silently transposable (spec §3.1).
///
/// Unconstrained at the type level: the coordinate invariant `1 <= start <= end` is checked
/// where a span is built from detector output, the same place [`crate::ng::region_typing::segment_criteria::SsrSegment`]
/// checks it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct TractSpan {
    /// First base of the tract, 1-based inclusive.
    pub start: Position,
    /// Last base of the tract, 1-based inclusive.
    pub end: Position,
}

impl TractSpan {
    /// The tract's length in bases. Inclusive on both ends, so a one-base span is 1.
    ///
    /// Assumes `start <= end`; a reversed span is a bug at the site that built it, and the
    /// `debug_assert` is where it is caught rather than in a length that silently wraps.
    #[inline]
    pub fn len_bp(self) -> u64 {
        debug_assert!(
            self.start.get() <= self.end.get(),
            "reversed span {}..={}",
            self.start.get(),
            self.end.get()
        );
        self.end.get() - self.start.get() + 1
    }

    /// How many whole motif copies of `period` bases the span holds.
    ///
    /// This is the **repeat count**: half of a stratum's key, and what a copy floor is
    /// compared against. It is derived rather than stored, so it cannot disagree with the
    /// span it comes from (spec §3.1).
    #[inline]
    pub fn repeat_count(self, period: u8) -> u64 {
        debug_assert!(period >= 1, "period 0 has no meaning");
        self.len_bp() / u64::from(period)
    }
}

/// One tandem repeat as the catalog holds it: what the scanner found, where the
/// whole-motif cut falls, and the two values that needed the bases.
///
/// **Both spans are stored because classification uses both, at different steps** (spec
/// §3.1): bundling and the pre-screen run on [`detected`](Self::detected), while the copy
/// floor, the purity floor, the flank test and the emitted locus all use
/// [`trimmed`](Self::trimmed). A row carrying one would send its reader back to the FASTA
/// for the other, which is the cost this file exists to remove.
#[derive(Debug, Clone, PartialEq)]
pub struct FoundRepeat {
    /// Which contig, as an index into the header's contig table.
    pub contig: ContigId,
    /// The span the scanner reported, converted from the detector's 0-based half-open
    /// [`crate::ng::tandem_repeat::RepeatInterval`] at the builder's edge.
    pub detected: TractSpan,
    /// The same tract cut back to whole motif copies at both ends, or `None` when no clean
    /// cut exists (`minimal_trim` fails).
    ///
    /// **This is the locus.** `None` is a real state the file records rather than a row it
    /// drops: such a tract can still be a **bundle member** for its neighbour, and dropping
    /// it would change how the neighbour is classified.
    pub trimmed: Option<TractSpan>,
    /// Motif length in bases — the axis a reader filters on, and half of a stratum's key.
    pub period: u8,
    /// The Ruzzo–Tompa segment total the detector scored this tract at.
    ///
    /// **Not a TRF score**, and not on the same scale: a threshold carried across from TRF
    /// would not mean here what it meant there
    /// (`crate::ng::region_typing::segment_criteria::SsrSegmentCriteria::min_score`).
    pub score: i32,
    /// The repeat unit, verbatim and phase-faithful — the tract's first `period` bases,
    /// uppercased.
    pub motif: Motif,
    /// Fraction of the **trimmed** tract matching a perfect motif tiling, in `[0, 1]`.
    /// `None` exactly when [`trimmed`](Self::trimmed) is.
    pub purity: Option<f32>,
}

impl FoundRepeat {
    /// The stratum this repeat belongs to — `(period, repeat count)` — or `None` when the
    /// tract has no clean trim and so is no locus under any policy.
    ///
    /// The count is measured on the **trimmed** span, which is the count `classify`
    /// applies its copy floor to (spec §3.1).
    #[inline]
    pub fn stratum(&self) -> Option<(u8, u64)> {
        self.trimmed
            .map(|span| (self.period, span.repeat_count(self.period)))
    }
}

/// What a catalog file says about itself, readable without touching a single row.
///
/// Stored in the file's footer metadata (spec §3.4, §3.5). Its job is to answer two
/// questions before any row is trusted: *is this the reference I am working on?* and *was
/// it built permissively enough to answer what I am about to ask?*
#[derive(Debug, Clone, PartialEq)]
pub struct RepeatCatalogHeader {
    /// One entry per contig, in reference order — name, length and the per-contig MD5 (the
    /// `@SQ M5`). The MD5s are what "is this the right reference?" is answered with, and
    /// they name **which** contig differs; the lengths are load-bearing too, because the
    /// flank test asks about the contig's own end.
    pub contigs: Vec<ContigInfo>,
    /// The whole-reference digest: every contig's uppercased bases concatenated in file
    /// order. Being an in-order concatenation, it also pins the contig order.
    pub reference_md5: [u8; 16],
    /// The criteria the builder was given. Only its period range, copy floors and minimum
    /// flank were **applied**; the rest are recorded for provenance, because a reader
    /// chooses those freely (spec §4.2).
    pub built_under: StrRepeatCriteria,
    /// The scoring weights every stored tract came out of. **Equal or refuse**: other
    /// weights are a different set of tracts, not a subset of this one (spec §4.2).
    pub scan: ScanParams,
    /// The crate version that wrote the file, so a change in the detector invalidates it
    /// even when every setting matches.
    pub tool_version: String,
    /// The longest tract stored, per contig, in the header's contig order — 0 where a
    /// contig holds no repeats.
    ///
    /// **This is what makes a windowed read exact.** Classification is not local: a repeat
    /// outside the stretch you asked for still bundles with one inside it, and a long array
    /// outside still swallows a locus inside. So a reader that wants a stretch must also
    /// read every row that could reach into it, and the reach of the longest row is the
    /// bound. Guessing it would be a silent wrong answer near long arrays; the builder
    /// knows it exactly and it costs one number per contig.
    pub longest_tract_bp: Vec<u64>,
}

impl RepeatCatalogHeader {
    /// How far outside a requested stretch a row can sit and still change the answer inside
    /// it: the contig's longest tract, plus the bundle radius the reader is asking with.
    ///
    /// A read window padded by this holds every row that can affect the answer, and no
    /// policy the file serves can need more — a wider bundle radius is the reader's to pass
    /// in, and a longer tract than the longest one stored does not exist in that contig.
    pub fn reach_into_window(&self, contig: ContigId, bundle_threshold: u64) -> u64 {
        let longest = self
            .longest_tract_bp
            .get(contig.get() as usize)
            .copied()
            .unwrap_or(0);
        longest + bundle_threshold
    }
}

/// Everything that can stop a catalog being built or read.
///
/// Every variant is fatal and every variant names what differed. None of them is a
/// warning, a fallback or a silent rebuild: a short answer here is a genome-wide count
/// taken over the wrong genome, and nothing downstream would notice (spec §4.3, §6).
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum RepeatCatalogError {
    /// No catalog at that path. **The one error a caller can act on**, so it names the
    /// command that writes one — distinct from a catalog that is present but wrong.
    #[error(
        "no repeat catalog at {path}; build one with `pop_var_caller_exp repeat-catalog --reference <ref.fa>`"
    )]
    NotFound {
        /// Where a catalog was looked for.
        path: PathBuf,
    },

    /// A contig's MD5 in the file differs from the one this run computed from the FASTA:
    /// the catalog describes a different reference.
    #[error(
        "catalog does not describe this reference: contig {contig} is {found} here, catalog says {expected}"
    )]
    ContigDigestMismatch {
        /// The contig whose bases differ.
        contig: String,
        /// The digest this run computed from the FASTA, hex.
        found: String,
        /// The digest the catalog's header carries, hex.
        expected: String,
    },

    /// The contig **tables** disagree — a contig present in one and not the other, a
    /// length that differs, or the same contigs in a different order. Order matters
    /// because rows address contigs by index.
    #[error("catalog contig table does not match the reference: {detail}")]
    ContigTableMismatch {
        /// Which contig, and how the two tables differ there.
        detail: String,
    },

    /// The reader asked about tracts the file does not hold — a copy floor below the built
    /// one, a flank below it, or a period outside its range (spec §4.1).
    #[error("catalog is not permissive enough: {0}")]
    CriteriaTooPermissive(#[from] CriteriaRefusal),

    /// The reader wants tracts scored with different weights. Not a filter over the file:
    /// the weights decide where each segment starts and stops, so a tract can grow, shrink
    /// or vanish (spec §4.2).
    #[error(
        "catalog was scored with match {built_match_reward}/mismatch {built_mismatch_penalty}, \
         reader asked for match {wanted_match_reward}/mismatch {wanted_mismatch_penalty}"
    )]
    ScoringWeightsDiffer {
        /// Match reward the file was built with.
        built_match_reward: i32,
        /// Mismatch penalty the file was built with.
        built_mismatch_penalty: i32,
        /// Match reward the reader asked for.
        wanted_match_reward: i32,
        /// Mismatch penalty the reader asked for.
        wanted_mismatch_penalty: i32,
    },

    /// The file was written by a different version of the detector, so its rows are not
    /// necessarily the ones this build would produce.
    #[error("catalog was written by tool version {built}, this is {running}")]
    ToolVersionDiffers {
        /// Version recorded in the file.
        built: String,
        /// Version running now.
        running: String,
    },

    /// The file cannot be opened, has no footer (a build that died mid-write), or is not a
    /// catalog at all.
    #[error("catalog file {path} is unreadable or truncated")]
    Unreadable {
        /// The file that could not be read.
        path: PathBuf,
        /// What the underlying reader said.
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },

    /// Writing the catalog failed.
    #[error("writing the catalog to {path}")]
    Write {
        /// The file being written.
        path: PathBuf,
        /// What the filesystem said.
        #[source]
        source: std::io::Error,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    fn span(start: u64, end: u64) -> TractSpan {
        TractSpan {
            start: Position(start),
            end: Position(end),
        }
    }

    #[test]
    fn a_span_is_inclusive_at_both_ends() {
        assert_eq!(
            span(10, 10).len_bp(),
            1,
            "a one-base tract is one base long"
        );
        assert_eq!(span(10, 19).len_bp(), 10);
    }

    #[test]
    fn the_repeat_count_is_whole_copies_only() {
        // 10 bases of a period-3 motif is three whole copies and a base left over.
        assert_eq!(span(1, 10).repeat_count(3), 3);
        assert_eq!(span(1, 12).repeat_count(3), 4);
        assert_eq!(span(1, 5).repeat_count(1), 5);
    }

    #[test]
    fn a_repeat_with_no_clean_trim_belongs_to_no_stratum() {
        let untrimmable = FoundRepeat {
            contig: ContigId(0),
            detected: span(100, 115),
            trimmed: None,
            period: 3,
            score: 20,
            motif: Motif::new(b"CAG").expect("valid motif"),
            purity: None,
        };
        assert_eq!(untrimmable.stratum(), None);
    }

    #[test]
    fn the_stratum_counts_the_trimmed_tract_not_the_detected_one() {
        // Detected 100..=115 (16 bases, 5 copies of a 3-mer); trimmed to 102..=113
        // (12 bases, 4 copies). The stratum must follow the trim, because that is the
        // count `classify` applies its copy floor to (spec §3.1).
        let repeat = FoundRepeat {
            contig: ContigId(0),
            detected: span(100, 115),
            trimmed: Some(span(102, 113)),
            period: 3,
            score: 20,
            motif: Motif::new(b"CAG").expect("valid motif"),
            purity: Some(1.0),
        };
        assert_eq!(repeat.stratum(), Some((3, 4)));
    }
}
