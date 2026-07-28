//! ng's generic locus generator — the pileup walk, copied from production.
//!
//! A **folder**, unlike its `ssr.rs` sibling, because it holds the whole copied
//! walker plus (later) the generator that wraps it. Production is not edited:
//! ng copies rather than reaches in, so there is no visibility lift and no field
//! on a frozen type (`doc/devel/ng/spec/locus_generation_pileup.md` §3,
//! `doc/devel/ng/arch/locus_generation_pileup.md` *Module home*).
//!
//! # Seven files, verbatim
//!
//! [`genome_walk`], [`open_record`], [`cigar_cursor`], [`decompose`],
//! [`active_read_set`], [`chain_id_allocator`] and [`errors`] are transcribed
//! from `src/pileup/walker/` **unchanged**, and are proven to compute exactly
//! what production computes before a line of them is edited (spec §3, the
//! stage-1 differential). The rule that paid three times on this branch is
//! *transcribe first, change second*: a copy that is provably production is the
//! baseline every later change is measured against, and without it the
//! generator's deliberate divergences could not be told from transcription
//! slips.
//!
//! So the copies still emit production's [`PileupRecord`](crate::pileup_record::PileupRecord),
//! not ng's `SampleLocusObservations`. That changes in the next plan.
//!
//! **One file is renamed on the way in: `driver.rs` → [`genome_walk`]** — it is
//! the only one of the seven named for a *role* rather than for what it owns,
//! and "driver" answers *driver of what?* with nothing. `genome_walk` names the
//! one job that file has and the others do not: advancing a position cursor
//! along genome coordinates over an active read set. The **type** keeps
//! production's `PileupWalker`, so the differential reads as a straight
//! comparison. *(A walk covers one region, not the genome; `genome_walk` names
//! the axis it advances along, not the extent.)*
//!
//! # What this module re-exports, and why
//!
//! The copied files reach their shared vocabulary through `super::` — that is
//! how they were written against `pileup/walker/mod.rs`, and leaving those paths
//! alone is what keeps the transcription verbatim. This module therefore stands
//! in for production's `walker/mod.rs`, drawing each name from wherever it now
//! lives:
//!
//! - [`PreparedRead`], [`MateRole`], [`ReadLengthError`] — **ng's**, copied and
//!   extended with `read_group` (spec §6). They are the reason the whole walker
//!   is copied: every one of the seven names `PreparedRead` in its signatures.
//! - [`CigarOp`], [`WalkerConfig`] and the `DEFAULT_*` constants — **production's,
//!   reused as-is**. ng does not modify them, so it does not copy them; the
//!   constants are reached by name rather than by literal so there is one source
//!   of truth until ng deliberately diverges.

mod active_read_set;
mod chain_id_allocator;
mod cigar_cursor;
mod decompose;
mod errors;
mod genome_walk;
mod open_record;

// The vocabulary the copied files resolve through `super::`. Production's own
// `walker/mod.rs` declares exactly these names around exactly these modules.
pub use crate::ng::read::prepared_read::{MateRole, PreparedRead, ReadLengthError};
pub use crate::pileup::walker::{
    CigarOp, DEFAULT_MATE_LOOKUP_WINDOW, DEFAULT_MAX_INDEL_COLUMN_DEPTH, DEFAULT_MAX_RECORD_SPAN,
    DEFAULT_MAX_SNP_COLUMN_DEPTH, WalkerConfig,
};

pub use chain_id_allocator::DEFAULT_MAX_ACTIVE_READS;
pub use errors::WalkerError;
pub use genome_walk::{PileupWalker, RunSummary, run};
