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
//! observations in parallel); and steps 6 to 9, the calling loop and what it drives
//! ([`calling`] — so far the vocabulary its four sub-modules will share, and two of the
//! four: step 8's [`calling::genotype_prior`], as its folder and the four files the plan
//! fills, and step 6's [`calling::allele_candidates`], so far the two constants its
//! admission rule is made of); and step 11a's hidden-duplication filter
//! ([`paralog`] — production's statistics copied in, all four of them: the per-sample
//! coverage model that says what one copy's depth looks like, the per-locus score that
//! weighs a collapsed pair of gene copies against a real variant, how common hidden
//! duplications are in this run, and the false-discovery curve that turns the operator's
//! target into a cut).
//!
//! **Production is frozen.** ng is a from-scratch caller: it does not edit
//! `src/ssr/` or `src/regions.rs` — nor, since the generic locus generator's port,
//! `src/pileup/`, `src/psp/`, `src/var_calling/` or `src/vcf/` — and it does not
//! depend on trf-mod. Where ng needs a production behaviour in a different shape, it
//! **copies the code and changes its own version** (owner, 2026-07-16); reuse is for
//! what costs production nothing. Winning steps are ported back only after the
//! experiments ng exists to run have decided something.
//!
//! **A test may read production as an oracle, and a handful do** — `scanner_parity`
//! against `src/ssr/` and `calling::genotype_table_parity` against `src/var_calling/`
//! were the first two; there are now several more, and the way to find them is
//! `grep -rnE 'use crate::|include_str!\("\.\./\.\.' src/ng | grep -v 'crate::ng'` rather
//! than a list here that goes stale. **The second alternative matters**: the copy guards
//! (`paralog/copy_fidelity.rs`, `locus_generation/pileup/copy_fidelity.rs`) read production's
//! source as *text* at compile time rather than importing from it, and a `use`-only sweep
//! does not see them at all. Every one is `#[cfg(test)]`, so nothing shipped depends on production; the
//! direction that matters is the other one, and production still depends on nothing in ng.
//! A port's whole claim is that it agrees with what it was ported from, and only production
//! can settle that. **Every occurrence outside a `#[cfg(test)]` module needs a stated
//! reason** — `paralog::coverage_model` has the one that exists today, and its own header
//! gives the reason and the date it ends.
//!
//! **One such oracle cost production one line, and it is the only edit ng has made to
//! a frozen tree.** `posterior_engine.rs` declared `mod shape;` privately, which put
//! `GenotypeShape` — the thing `calling::genotype_table` is a port *of* — out of reach
//! of any test that could check the port. It is now `pub(crate)` (owner, 2026-08-21).
//! No behaviour moved and nothing was re-exported. The rule this bends is worth stating
//! precisely rather than quietly: **ng may widen a production item's visibility so a
//! parity test can see it, and may change nothing else**; anything that would alter what
//! production computes is still a copy-into-ng, not an edit.
//!
//! **The mirror rule, and it is a different one with a different subject: ng may widen its
//! own copy's visibility where production's module was private and ng's is not.**
//! [`paralog::calibration`] is a span of `src/var_calling/paralog_filter/calibrate.rs`, which
//! production keeps `pub(crate)` inside a private module; ng re-exports the same names from
//! [`paralog`], and a `pub use` cannot re-export a `pub(crate)` item, so five lines read
//! `pub` where production reads `pub(crate)`. **Production is untouched.** Each such line is
//! declared and checked in `paralog/copy_fidelity.rs`, and the check is exact: the two lines
//! must be identical once `pub(crate)` becomes `pub`.
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
pub mod calling;
pub mod genetics;
pub mod locus_generation;
pub mod paralog;
pub mod parameter_estimation;
pub mod psp;
pub mod raw_chrom_reader;
pub mod read;
pub mod ref_seq;
pub mod reference_info;
pub mod region_typing;
pub mod repeat_catalog;
pub mod run;
pub mod segmentation_inputs;
pub mod tandem_repeat;
pub mod types;
pub mod vcf;
pub mod window_coverage;

pub use ref_seq::{
    ContigTable, EvictableRefSeq, InMemoryRefSeq, RawRefSeq, RefSeq, RefSeqError, ResidentRefSeq,
    WindowedRefSeq,
};
pub use types::{
    BaseQual, Bp, ContigId, DomainError, ErrorRate, GenotypeFrequency, InbreedingF, MapQual,
    MismatchFraction, Ploidy,
};
