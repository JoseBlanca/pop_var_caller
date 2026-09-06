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
pub use spill::{SpillEntry, SpillError, SpillReader, SpillWriter, SpilledSample};
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

/// One sample's coverage at one locus: the GC fraction of the window centred on it and that
/// window's mean read depth. **Both fields `NaN` where the sample has no usable window
/// there** — stored and compared by bit pattern, never by `==`.
///
/// **This is a stand-in, and it is scheduled for deletion.** The type belongs to the window
/// coverage design, which puts it in `src/ng/window_coverage/` on branch `ng-window-coverage`;
/// neither that module nor its specification is on `main`, and this module must not create or
/// touch anything under that path. So the two fields that design pins are declared here, and
/// this declaration goes away when the branches meet — the codec reads the two fields and
/// nothing else, so the swap is an import.
///
/// The crate holds a second, unrelated `WindowCoverage` in
/// [`crate::sample_summary::coverage`], which is production's per-tile summary; the two do not
/// meet, and this one outlives the rebase in neither name nor place.
#[derive(Clone, Copy, Debug)]
pub struct WindowCoverage {
    /// The share of the window's bases that are G or C, between 0 and 1 — or `NaN` where the
    /// sample has no usable window at this locus.
    pub gc_fraction: f32,
    /// The window's mean read depth in this sample, or `NaN` in the same case.
    pub mean_depth: f32,
}
