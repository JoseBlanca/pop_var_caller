//! Read-level pipeline stages — turning a decoded alignment record into locus
//! evidence. Two steps share this module because they share the same reference
//! accessor and read filtering's output (`MappedRead`) is read preparation's
//! input:
//!
//! - **step 1 — [`filtering`]** — the fixed whole-read keep/drop prelude. A
//!   single file, no bake-off: it is a fixed prelude with no competing
//!   implementations.
//! - **step 2 — read preparation ([`ReadPreparer`])** — canonicalise the line-up
//!   the mapper gave a kept read, against the reference around that read's own
//!   span. **Generic path only:** the STR path re-aligns every spanning read
//!   against its tract, so canonicalizing the mapper's CIGAR there would be work
//!   nothing reads (`doc/devel/ng/spec/read_preparation.md` §1). Its swappable
//!   implementations sit here as siblings — [`left_align`] today, a
//!   pass-through-only null arm and the gated re-align mode later.
//!
//! This is the deliberate, documented deviation from the "one folder per step"
//! rule in `doc/devel/ng/arch/module_layout.md` principle 1 — see
//! `doc/devel/ng/spec/read_filtering.md` §2.1.
//!
//! [`input`] joins them as step 1's **input edge** rather than a step of its
//! own: opening and validating an alignment file, serving a region as an
//! ordered stream, and merging a sample's several files into one. It feeds
//! [`filtering`] — its region readers are the `RecordSource` the filter
//! consumes — so it belongs beside it (`module_layout.md` principle 1, note b).

pub mod filtering;
pub mod input;
pub mod left_align;

pub use filtering::{
    BamRecordSource, CramRecordSource, NoodlesRawRecord, RawRecord, ReadFilter, ReadFilterConfig,
    ReadFilterCounts, ReadFilterError, RecordSource,
};

use crate::bam::alignment_input::MappedRead;
use crate::ng::ref_seq::RefSeqError;
use crate::pileup::walker::PreparedRead;

/// A fatal, run-level failure from read preparation. **Never a per-read verdict:** a read
/// that yields no usable observation is `Ok(None)`, and the run continues (spec §7).
///
/// In v1 there is exactly one source — the reference fetch that left-alignment needs, which
/// only indel-bearing reads reach at all. A failed fetch means a malformed record or a
/// mis-built reference, so it ends the run rather than silently leaving one read
/// un-normalised, which is what production does
/// ([read_processor.rs](../../../pileup/per_sample/read_processor.rs)) and what ng
/// deliberately diverges from (spec §7, §9). Step 1 made the same call for the same reason
/// ([`ReadFilterError::Reference`](filtering::ReadFilterError)).
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum ReadPrepError {
    /// The reference window an indel-bearing read needed could not be fetched.
    #[error("reference access failed during read preparation")]
    Reference(#[source] RefSeqError),
}

/// Canonicalise the line-up a mapper gave one read, against the reference around that read's
/// **own** span — the generic path's step 2.
///
/// # What every implementation upholds
///
/// **Per read, and independent of every other read.** No mate, no neighbouring read, no
/// locus-column context. That is what makes preparation parallel with deterministic output,
/// and why the genuinely pairwise work — reconciling overlapping mates' qualities — sits
/// downstream in the pileup walker. It is also why local reassembly cannot be an
/// implementation here: assembling haplotypes needs every read in a region at once.
///
/// **Locus-independent**, so one output serves *every* locus the read overlaps; hence no
/// locus argument.
///
/// **It re-places what the read already says; it never invents sequence.**
///
/// # The two kinds of bad news, which must not be confused
///
/// `Ok(None)` is a read that produced no usable observation — counted by reason, run
/// continues. `Err` is a broken run. An `Option` alone has one channel and could express only
/// one of them, which is why the return type carries both (spec §7). **No implementation in
/// v1 returns `Ok(None)`:** the only step that could decline a read was BAQ, deferred sine die.
/// The variant is kept for the implementations that will reintroduce a decline.
///
/// # Ownership and dispatch
///
/// The read arrives **by value** so its `seq`, `qual` and `cigar` buffers *move* into the
/// [`PreparedRead`] instead of being cloned on a path that runs for every read of every
/// sample; step 1 already yields owned reads, so this costs its caller nothing.
/// Implementations are selected by **generic type parameter, never behind `dyn`** — a virtual
/// call per read is a cost this path cannot carry.
///
/// **No reference argument:** an implementation that needs a reference holds its own accessor
/// as a field, and one that does not, does not — so "this preparer needs no reference" is a
/// compile-time fact rather than a runtime `Option` with an undefined case (spec §6).
///
/// Design: `doc/devel/ng/spec/read_preparation.md` §6, `doc/devel/ng/arch/read_preparation.md` §2.
pub trait ReadPreparer {
    /// Reused buffers — allocated once per worker, never per read. `()` for an implementation
    /// that needs none.
    ///
    /// `Default` is the only way a caller builds one, so whatever it picks must be invisible
    /// at every call site: scratch is **buffers only**, empty and grown on demand, and inert
    /// with respect to the output. An implementation needing a parameter that changes a result
    /// puts it in its own constructor state, where it is visible.
    type Scratch: Default;

    /// Prepare one filtered read.
    ///
    /// `Ok(Some(read))` is the prepared read; `Ok(None)` is "no usable observation", tallied
    /// by the caller; `Err` ends the run. An implementation that allocates per call — rather
    /// than reusing `scratch` — is a defect, not a slow path.
    fn prepare_read(
        &self,
        read: MappedRead,
        scratch: &mut Self::Scratch,
    ) -> Result<Option<PreparedRead>, ReadPrepError>;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::types::ContigId;
    use crate::pileup::per_sample::baq_engine::prepare_passthrough;
    use crate::pileup::walker::CigarOp;

    /// A minimal mapped read: 4 matched bases at position 10 on contig 0.
    fn mapped_read() -> MappedRead {
        MappedRead {
            qname: b"read1".to_vec(),
            flag: 0,
            ref_id: 0,
            pos: 10,
            mapq: 60,
            cigar: vec![CigarOp::Match(4)],
            seq: b"ACGT".to_vec(),
            qual: vec![30, 30, 30, 30],
            mate_ref_id: None,
            mate_pos: None,
            adaptor_boundary: None,
            source_file_index: 0,
        }
    }

    /// Prepares every read, by handing the whole build to production's `--no-baq` path.
    struct AlwaysPrepares;
    impl ReadPreparer for AlwaysPrepares {
        type Scratch = ();
        fn prepare_read(
            &self,
            read: MappedRead,
            _scratch: &mut (),
        ) -> Result<Option<PreparedRead>, ReadPrepError> {
            Ok(Some(prepare_passthrough(read, 0)))
        }
    }

    /// Declines every read — the `Ok(None)` channel no v1 implementation uses, kept
    /// implementable so the shape is pinned before BAQ or the re-align mode need it.
    struct AlwaysDeclines;
    impl ReadPreparer for AlwaysDeclines {
        type Scratch = ();
        fn prepare_read(
            &self,
            _read: MappedRead,
            _scratch: &mut (),
        ) -> Result<Option<PreparedRead>, ReadPrepError> {
            Ok(None)
        }
    }

    /// Fails the run — the channel that exists because a decline and a broken run must not
    /// share one.
    struct AlwaysFails;
    impl ReadPreparer for AlwaysFails {
        type Scratch = ();
        fn prepare_read(
            &self,
            _read: MappedRead,
            _scratch: &mut (),
        ) -> Result<Option<PreparedRead>, ReadPrepError> {
            Err(ReadPrepError::Reference(RefSeqError::UnknownContig(
                ContigId(7),
            )))
        }
    }

    /// A scratch that is **not** `()`, so `Scratch: Default` has something to prove.
    #[derive(Default)]
    struct CountingScratch {
        buffer: Vec<u8>,
        calls: usize,
    }

    /// Reuses its scratch across reads rather than allocating per call — the contract the
    /// associated type exists for.
    struct ReusesScratch;
    impl ReadPreparer for ReusesScratch {
        type Scratch = CountingScratch;
        fn prepare_read(
            &self,
            read: MappedRead,
            scratch: &mut CountingScratch,
        ) -> Result<Option<PreparedRead>, ReadPrepError> {
            scratch.calls += 1;
            scratch.buffer.clear();
            scratch.buffer.extend_from_slice(&read.seq);
            Ok(Some(prepare_passthrough(read, 0)))
        }
    }

    /// Drive a preparer **through a generic bound**, which is the only place `Scratch: Default`
    /// and the never-`dyn` rule actually bite: a concrete call compiles whether or not the trait
    /// requires them, so a test using concrete types alone would pin neither. Here the scratch
    /// comes from `Default` and is threaded across two reads, so both the bound and the reuse
    /// contract are load-bearing — drop either and this helper stops compiling.
    fn prepare_pair<P: ReadPreparer>(
        preparer: &P,
        reads: [MappedRead; 2],
    ) -> Result<Vec<Option<PreparedRead>>, ReadPrepError> {
        let mut scratch = P::Scratch::default();
        reads
            .into_iter()
            .map(|read| preparer.prepare_read(read, &mut scratch))
            .collect()
    }

    #[test]
    fn a_preparer_yields_the_prepared_read_through_a_generic_caller() {
        let prepared = prepare_pair(&AlwaysPrepares, [mapped_read(), mapped_read()])
            .expect("this preparer never fails");
        assert_eq!(prepared.len(), 2);
        let first = prepared[0].as_ref().expect("prepared");
        assert_eq!(first.alignment_start, 10);
        assert_eq!(first.seq, b"ACGT");
    }

    /// The two `Ok` outcomes are distinguishable, and neither is an error. `Ok(None)` means the
    /// run continues without this read; only `Err` ends it.
    #[test]
    fn declining_a_read_is_an_answer_not_an_error() {
        let outcomes =
            prepare_pair(&AlwaysDeclines, [mapped_read(), mapped_read()]).expect("not an error");
        assert!(outcomes.iter().all(Option::is_none));
    }

    /// A broken run surfaces as `Err`, never folded into `Ok(None)` — the distinction §7 exists
    /// to enforce.
    #[test]
    fn a_broken_run_is_an_error_not_a_decline() {
        let outcome = prepare_pair(&AlwaysFails, [mapped_read(), mapped_read()]);
        assert!(matches!(
            outcome,
            Err(ReadPrepError::Reference(RefSeqError::UnknownContig(
                ContigId(7)
            )))
        ));
    }

    /// One scratch serves every read: it is built once by the generic caller and threaded
    /// through, so an implementation can reuse its buffers instead of allocating per read.
    #[test]
    fn one_scratch_is_reused_across_reads() {
        let preparer = ReusesScratch;
        let mut scratch = CountingScratch::default();
        assert_eq!(scratch.calls, 0);
        for _ in 0..3 {
            preparer
                .prepare_read(mapped_read(), &mut scratch)
                .expect("never fails");
        }
        assert_eq!(scratch.calls, 3, "the same scratch saw every read");
        assert_eq!(scratch.buffer, b"ACGT", "and it carries state between them");
    }
}
