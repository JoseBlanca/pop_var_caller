//! Errors raised by the pileup walker. Each variant carries enough
//! context (qname, chromosome, position) to point at the offending
//! input from the error alone, per `ia/specs/design_principles.md`
//! principle 6 ("typed errors at module boundaries").
//!
//! **No longer a verbatim copy — A0 (plan 3).** Copied from
//! `src/pileup/walker/errors.rs`, then changed in exactly one place:
//! [`WalkerError::Fasta`]'s source is ng's [`RefSeqError`], because ng's walker
//! fetches through [`RefSeq`](crate::ng::ref_seq::RefSeq) rather than through
//! production's `MultiChromRefFetcher`. `copy_fidelity.rs` released this file in
//! that commit; everything else here is still production's, line for line.

use thiserror::Error;

use crate::ng::ref_seq::RefSeqError;

#[derive(Error, Debug)]
#[non_exhaustive]
pub enum WalkerError {
    #[error(
        "out-of-order read: qname='{qname}' at (chrom_id={chrom_id}, pos={pos}) regresses from \
         (chrom_id={prev_chrom_id}, pos={prev_pos})"
    )]
    OutOfOrder {
        qname: String,
        prev_chrom_id: u32,
        prev_pos: u32,
        chrom_id: u32,
        pos: u32,
    },

    #[error(
        "read decoded zero reference bases: qname='{qname}' at (chrom_id={chrom_id}, pos={pos})"
    )]
    ZeroRefSpan {
        qname: String,
        chrom_id: u32,
        pos: u32,
    },

    /// **Unreachable from ng's walk, and kept because the type is production's.**
    /// `WalkerState::new` hands the chain-id allocator a cap of `u32::MAX`; the real
    /// ceiling, `max_active_reads`, is enforced in `admit_read` by refusing the read.
    /// Production's walker, which shares this variant's text through its own copy of
    /// `errors.rs`, still raises it.
    #[error(
        "active-read cap exceeded (cap={cap}) at chrom_id={chrom_id} pos={pos}; \
         consider raising --max-active-reads or pre-filtering this region"
    )]
    ActiveReadsExhausted { cap: u32, chrom_id: u32, pos: u32 },

    #[error(
        "phase-chain id space exhausted at chrom_id={chrom_id} pos={pos}: \
         this .psp file has reached 2^64 unique read identifiers, the per-file \
         limit imposed by the u64 chain id encoding"
    )]
    ChainIdSpaceExhausted { chrom_id: u32, pos: u32 },

    #[error(
        "pending-mates map exceeded its defensive cap (cap={cap}) at \
         chrom_id={chrom_id} pos={pos}: more than {cap} reads are waiting for a mate \
         at once. Two inputs do this. The common one is depth — every read whose \
         partner lies further along waits here, so a region deep enough fills the map \
         with correctly-paired reads; lower --max-active-reads so the walk holds fewer \
         reads open. The other is broken pairing (every paired read flagged \
         FirstOfPair, no SecondOfPair ever arriving), which shows up as a \
         mate_lookup_evictions count near the read count rather than near zero."
    )]
    PendingMatesExhausted { cap: usize, chrom_id: u32, pos: u32 },

    #[error(
        "open record reference span exceeded MAX_RECORD_SPAN: anchor (chrom_id={chrom_id}, \
         pos={pos}) reached span {span} (cap={cap}); upstream filter should have rejected the \
         underlying read"
    )]
    RecordTooWide {
        chrom_id: u32,
        pos: u32,
        span: u32,
        cap: u32,
    },

    #[error("FASTA fetch failed at chrom_id={chrom_id} for [{start}, {start_plus_len}): {source}")]
    Fasta {
        chrom_id: u32,
        start: u32,
        start_plus_len: u32,
        #[source]
        source: RefSeqError,
    },

    #[error(
        "internal invariant violated: {detail} (qname='{qname}' chrom_id={chrom_id} pos={pos})"
    )]
    Internal {
        detail: String,
        qname: String,
        chrom_id: u32,
        pos: u32,
    },

    #[error("malformed PreparedRead at qname='{qname}' (chrom_id={chrom_id}, pos={pos}): {reason}")]
    MalformedRead {
        reason: String,
        qname: String,
        chrom_id: u32,
        pos: u32,
    },
}
