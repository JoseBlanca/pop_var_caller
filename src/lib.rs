//! # `pop_var_caller`
//!
//! A population variant caller: SNPs, indels and repeat tracts, from one sample to a cohort, from
//! aligned reads to one VCF. Its design is under `doc/devel/ng/` (`spec/` says what each part must
//! do, `arch/` how it is built).
//!
//! The modules, roughly in the order a run reaches them:
//!
//! - [`bam`], [`fasta`], [`regions`] — the inputs: alignment files, the reference and region lists.
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
//! [`error_render`] and [`iter_ext`] are small helpers the rest shares. The command line is
//! [`pop_var_caller_exp`].
//!
//! The tandem-repeat scanner is this crate's own; nothing here depends on `trf-mod`, whose output
//! survives only as a committed test catalog (the test-only `golden_catalog`).
//!
//! **This caller was written as "ng", beside an older one, and replaced it on 2026-09-13.** The
//! older caller — "production" in this tree's comments — was a two-phase engine: a pileup, a
//! `.psp` store in a format this crate's [`psp`] does not share, and a cohort caller, under what
//! were `src/pileup/`, `src/psp/`, `src/var_calling/`, `src/ssr/`, `src/vcf/` and their neighbours
//! (`src/psp/` and `src/vcf/` now hold this crate's own modules of those names). ng
//! copied production's code wherever it needed one of its behaviours, changed its own copy, and
//! checked the copies against production in tests. Promotion Milestone C froze every answer those
//! tests compared against into fixtures — production's answers as its code stood at commit
//! `d9e7b076` — and Milestone D deleted production
//! (`doc/devel/implementation_plans/promote_ng_to_production.md`). So a comment that says
//! "production's", "ported from" or "copied from" means that caller, and its source is whole at
//! commit `74e424ea`, the last before the deletion. Milestone E moved ng's modules from `src/ng/`
//! up to the crate root; the design documents keep the name (`doc/devel/ng/`), and a comment that
//! says "ng" means this caller.
//!
//! ## Feature flags
//!
//! - `dhat-heap` — opt-in `dhat::Alloc` global allocator for heap
//!   profiling under benches and examples. Bench/example use only;
//!   not for production builds.
//! - `alloc-mimalloc` — the `mimalloc` global allocator, **on by
//!   default**: faster and smaller than the system allocator on this
//!   crate's workloads, measured on the cohort merge (and before
//!   that on the deleted production caller's `var-calling` path).
//!   `--no-default-features` opts out.
//!   Cannot hold the `#[global_allocator]` slot alongside `dhat-heap`,
//!   which wins it — so a heap profile is
//!   `--no-default-features --features dhat-heap`.

#![forbid(unsafe_code)]

#[cfg(test)]
pub(crate) mod golden_catalog;

#[cfg(test)]
pub(crate) mod production_post_filter;

#[cfg(test)]
mod scanner_parity;

pub mod alignment;
pub mod bam;
pub mod calling;
pub mod error_render;
pub mod fasta;
pub mod genetics;
pub mod iter_ext;
pub mod locus_generation;
pub mod paralog;
pub mod parameter_estimation;
pub mod pop_var_caller_exp;
pub mod psp;
pub mod raw_chrom_reader;
pub mod read;
pub mod ref_seq;
pub mod reference_info;
pub mod region_typing;
pub mod regions;
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
