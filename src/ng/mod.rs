//! ng — the caller: SNPs, indels and repeat tracts, from aligned reads to one VCF. Its design is
//! under `doc/devel/ng/` (`spec/` says what each part must do, `arch/` how it is built).
//!
//! The modules, roughly in the order a run reaches them:
//!
//! - [`reference_info`], [`ref_seq`], [`raw_chrom_reader`] — the reference: its contig table and
//!   digest, and its bases.
//! - [`tandem_repeat`], [`repeat_catalog`] — every tandem repeat in the reference, found once and
//!   written beside the FASTA.
//! - [`region_typing`], [`segmentation_inputs`] — what each stretch of the reference is — a *typed
//!   region*: a repeat tract or generic sequence — and the settings that decide it.
//! - [`read`], [`alignment`] — a decoded read turned into evidence: filtered, prepared,
//!   left-aligned, and lined up against the reference across repeat tracts.
//! - [`locus_generation`] — one sample's loci over those typed regions; [`psp`] — the per-sample
//!   store (one file per sample) the loci can be written to and the cohort merge, below, reads
//!   back.
//! - [`window_coverage`], [`paralog`] — each sample's depth and GC fraction by window, and the
//!   hidden-duplication filter's statistics built on them.
//! - [`parameter_estimation`] — the error rates, heterozygosity and inbreeding the caller runs on,
//!   measured from the cohort's own loci before anything is called.
//! - [`calling`], [`genetics`] — candidate alleles, read likelihoods, genotype priors, the loop
//!   that fits each site's allele frequencies, and site quality.
//! - [`run`] — what a run is made of: the walk over each sample's alignments, the psp files it can
//!   store, the parameter fit, the cohort merge that joins every sample's loci, and the calling
//!   that follows; [`vcf`] — the one file a run writes.
//! - [`types`] — the vocabulary the steps share.
//!
//! The tandem-repeat scanner is ng's own; nothing here depends on `trf-mod`, whose output survives
//! only as a committed test catalog (`golden_catalog`, test-only).
//!
//! The command line is [`crate::pop_var_caller_exp`].
//!
//! **ng began beside an older caller, and replaced it on 2026-09-13.** That caller —
//! "production" in this tree's comments — was a two-phase engine — a pileup, a `.psp` store in a
//! format ng's [`psp`] does not share, and a cohort caller — under `src/pileup/`, `src/psp/`,
//! `src/var_calling/`, `src/ssr/`, `src/vcf/` and their neighbours. ng copied its code wherever
//! it needed one of its behaviours and then changed its own copy, and it checked the copies against
//! production in tests. Promotion Milestone D
//! (`doc/devel/implementation_plans/promote_ng_to_production.md`) deleted production, after
//! Milestone C had frozen every answer those tests compared against into fixtures — production's
//! answers as its code stood at commit `d9e7b076`. So a comment here that says "production's", "ported from" or "copied from"
//! means that caller; its source is whole at commit `74e424ea`, the last before the deletion.

#[cfg(test)]
pub(crate) mod golden_catalog;

#[cfg(test)]
pub(crate) mod production_post_filter;

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
