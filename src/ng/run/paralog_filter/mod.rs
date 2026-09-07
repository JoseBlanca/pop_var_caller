//! **Where the hidden-duplication filter meets ng's run** — the file it parks called records
//! on, and the passes that fill and drain it.
//!
//! The filter's statistics live in [`crate::ng::paralog`] and are production's, copied. This
//! module is ng's own, and it exists because ng learns the filter's inputs later than
//! production does. Production fits each sample's coverage model when it opens its reader and
//! scores every record as it is called, so it reads its spill once. ng builds the coverage
//! histograms *during* the calling pass, so the models do not exist until that pass has
//! ended — which makes it two reads of the spill instead of one, and three passes over the
//! called records:
//!
//! 1. call, and write each finished record's line to a spill file beside the output;
//! 2. fit the models, score every spilled record, estimate how common duplications are in
//!    this run, and resolve the operator's target false-discovery rate to a cut;
//! 3. read the spill again and write the VCF, dropping or tagging the records above the cut.
//!
//! Spec: `doc/devel/ng/spec/hidden_paralog_filter.md`. Plan:
//! `doc/devel/ng/impl_plan/hidden_paralog_filter.md`.
//!
//! **What exists so far:** [`spill`]'s entry and its codec — what pass one writes and passes
//! two and three read back; [`spill_file`]'s [`SpillFile`], which says where those bytes live
//! and makes the file go away when the run ends, whatever way it ends; and [`patch`]'s
//! [`rewrite_filter_and_info`], which is how pass three puts the verdict on a line without
//! disturbing the columns it does not touch. Nothing fills a spill yet: the sink that appends to
//! it, the scoring context and the three passes are later steps of the plan above.

use crate::ng::types::GenomePosition;
use crate::ng::vcf::RecordPlace;

pub mod patch;
pub mod spill;
pub mod spill_file;

pub use patch::{LinePatchError, rewrite_filter_and_info};
pub use spill::{
    OnePositionSample, SpillEntry, SpillError, SpillReader, SpillWriter, SpilledSamples, WideSample,
};
pub use spill_file::{SpillFile, SpillFileError};

/// **Where pass three gets the place it hands the writer.**
///
/// [`VcfWriter::write_line`](crate::ng::vcf::VcfWriter::write_line) checks the order against the
/// place and never against the line's bytes, so a place typed out beside a line it does not
/// describe writes a VCF whose `POS` column runs backwards, with the check passing. Deriving it
/// from the entry the line came from is what makes that unbuildable: the three head fields are
/// read, not retyped.
impl From<&SpillEntry> for RecordPlace {
    fn from(entry: &SpillEntry) -> Self {
        Self {
            at: GenomePosition {
                contig: entry.contig,
                position: entry.position,
            },
            is_repeat_tract: entry.is_repeat_tract,
        }
    }
}

/// One sample's coverage at one locus, re-exported from the module that owns it.
///
/// **This was a stand-in until `ng-window-coverage` merged.** The two fields it declared —
/// the window's GC fraction and its mean read depth, both `NaN` where the sample has no usable
/// window — are exactly what [`crate::ng::window_coverage::WindowCoverage`] carries, so the
/// swap was this re-export and the deletion of the copy. The codec in [`spill`] reads the two
/// fields and nothing else, and did not change.
///
/// The crate holds a second, unrelated `WindowCoverage` in
/// [`crate::sample_summary::coverage`], which is production's per-tile summary; the two do not
/// meet.
pub use crate::ng::window_coverage::WindowCoverage;
