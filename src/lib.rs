//! # `pop_var_caller`
//!
//! Multi-sample variant caller — per-sample → cohort merge pipeline.
//! See `ia/specs/calling_pipeline_architecture.md` for the stage
//! breakdown.
//!
//! ## Feature flags
//!
//! - `dhat-heap` — opt-in `dhat::Alloc` global allocator for heap
//!   profiling under benches and examples. Bench/example use only;
//!   not for production builds.
//! - `alloc-mimalloc` — the `mimalloc` global allocator, **on by
//!   default**: faster and smaller than the system allocator on this
//!   crate's workloads, measured on both the production `var-calling`
//!   path and the ng cohort merge. `--no-default-features` opts out.
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
