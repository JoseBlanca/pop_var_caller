//! # `pop_var_caller`
//!
//! A population variant caller: SNPs, indels and repeat tracts, from one sample to a cohort, from
//! aligned reads to one VCF. The caller is [`ng`]; its command line is [`pop_var_caller_exp`].
//! [`bam`], [`fasta`] and [`regions`] read its inputs — alignment files, the reference and region
//! lists. They were first written for an older caller, deleted on 2026-09-13; [`ng`]'s header says
//! what that means for the comments in this tree.
//!
//! ## Feature flags
//!
//! - `dhat-heap` — opt-in `dhat::Alloc` global allocator for heap
//!   profiling under benches and examples. Bench/example use only;
//!   not for production builds.
//! - `alloc-mimalloc` — the `mimalloc` global allocator, **on by
//!   default**: faster and smaller than the system allocator on this
//!   crate's workloads, measured on the ng cohort merge (and before
//!   that on the deleted production caller's `var-calling` path). `--no-default-features` opts out.
//!   Cannot hold the `#[global_allocator]` slot alongside `dhat-heap`,
//!   which wins it — so a heap profile is
//!   `--no-default-features --features dhat-heap`.

#![forbid(unsafe_code)]

pub mod bam;
pub mod error_render;
pub mod fasta;
pub mod iter_ext;
pub mod ng;
pub mod pop_var_caller_exp;
pub mod regions;
