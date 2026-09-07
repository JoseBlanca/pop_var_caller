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
//! **What exists so far**, module by module:
//!
//! - [`spill`] — the entry a called record is parked as, and the codec that carries it.
//! - [`spill_file`] — [`SpillFile`]: where those bytes live, the three-stage life that stops a
//!   run writing over its own spill, and the `Drop` that unlinks it on every exit path.
//! - [`patch`] — [`rewrite_filter_and_info`]: how pass three puts the verdict on a line without
//!   disturbing the columns it does not touch.
//! - [`scoring_context`] — [`ParalogScoringContext`]: one fitted coverage model a sample, and the
//!   rule that turns a spilled row into the four numbers the scorer takes.
//! - [`pass_one`] — [`CalledRecordSink`]: the VCF writer while the filter is off, the spill while
//!   it is on.
//! - [`pass_two`] — [`score_the_parked_records_and_resolve_the_cut`]: read the spill, score every
//!   record, fit how common duplications are in this run, resolve the target to a cut.
//! - [`pass_three`] — [`write_the_records_the_filter_kept`]: read the spill again in step with
//!   those ratios and write the VCF, dropping or tagging what the cut removes.

use crate::ng::types::GenomePosition;
use crate::ng::vcf::RecordPlace;

pub mod finish;
pub mod pass_one;
pub mod pass_three;
pub mod pass_two;
pub mod patch;
pub mod scoring_context;
pub mod spill;
pub mod spill_file;

pub use finish::{
    FilteredRun, ParalogFilterError, WhatTheOperatorAskedFor, fit_score_and_write_the_calls,
    what_to_tell_the_operator,
};
pub use pass_one::{CalledRecordSink, PassOneError, SpillingSink, entry_for};
pub use pass_three::{
    HIDDEN_PARALOG_FILTER_ID, PassThreeError, WhatTheFilterDid, write_the_records_the_filter_kept,
};
pub use pass_two::{
    LrHistogramShape, NotATargetFdr, ParalogVerdicts, PassTwoError, TargetFdr,
    score_the_parked_records_and_resolve_the_cut,
};
pub use patch::{LinePatchError, rewrite_filter_and_info};
pub use scoring_context::{
    CohortSizeMismatch, CoverageFitConfigRefused, ParalogScoringContext, WhyNoCoverageModel,
};
pub use spill::{
    GenericLocusSample, RepeatTractSample, SpillEntry, SpillError, SpillReader, SpillWriter,
    SpilledSamples,
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
