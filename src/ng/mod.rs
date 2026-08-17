//! ng — the next-generation, step-decomposed variant/STR caller: a single-phase,
//! in-memory research lab (see `doc/devel/ng/`). Winning steps are ported back into the
//! production two-phase engine. Landed so far: the shared type vocabulary
//! ([`types`]); the `RefSeq` reference-sequence accessor ([`ref_seq`]); the step-1
//! read-filtering module ([`read`]); the tandem-repeat scanner primitive
//! ([`tandem_repeat`], a shared sequence primitive — types-and-scaffold stage);
//! step 3's typed-region generator ([`region_typing`] — Milestone A, the
//! segment-criteria port); the read-alignment module ([`alignment`] — a
//! folder of competing aligners, not a pipeline step); and step 4's locus
//! generation ([`locus_generation`] — the shared locus type and contract, the STR
//! generator, and the generic one's copied pileup walk); and step 4's parameter
//! pre-pass ([`parameter_estimation`] — the SNP/indel path, built milestone by
//! milestone); and the calling run ([`run`] — so far the cohort merge's parameters,
//! the first piece of the stage that turns the samples' observations into cohort
//! observations in parallel).
//!
//! **Production is frozen.** ng is a from-scratch caller: it does not edit
//! `src/ssr/` or `src/regions.rs` — nor, since the generic locus generator's port,
//! `src/pileup/`, `src/psp/`, `src/var_calling/` or `src/vcf/` — and it does not
//! depend on trf-mod. Where ng needs a production behaviour in a different shape, it
//! **copies the code and changes its own version** (owner, 2026-07-16); reuse is for
//! what costs production nothing. Winning steps are ported back only after the
//! experiments ng exists to run have decided something.
//!
//! **The heaviest instance of that rule so far is
//! [`locus_generation::pileup`]** — begun as a verbatim copy of `src/pileup/walker/`
//! (~5,500 lines), kept *provably* identical to its source so that ng's later,
//! deliberate divergences could be told apart from transcription slips. It is checked
//! textually, not asserted: `pileup/copy_fidelity.rs`.
//!
//! **Those divergences have now begun**, so "a verbatim copy" is no longer true of the
//! directory as a whole: three files have been released from the guard and changed on
//! purpose — the walker no longer fabricates the reference bases a read did not witness
//! — and four remain verbatim and guarded. `copy_fidelity.rs` names which are which, and
//! is the only thing that stays true as the balance shifts.

#[cfg(test)]
mod scanner_parity;

pub mod alignment;
pub mod locus_generation;
pub mod parameter_estimation;
pub mod raw_chrom_reader;
pub mod read;
pub mod ref_seq;
pub mod reference_info;
pub mod region_typing;
pub mod repeat_catalog;
pub mod run;
pub mod tandem_repeat;
pub mod types;

pub use ref_seq::{
    ContigTable, EvictableRefSeq, InMemoryRefSeq, RawRefSeq, RefSeq, RefSeqError, ResidentRefSeq,
    WindowedRefSeq,
};
pub use types::{
    BaseQual, Bp, ContigId, DomainError, ErrorRate, GenotypeFrequency, InbreedingF, MapQual,
    MismatchFraction, Ploidy,
};
