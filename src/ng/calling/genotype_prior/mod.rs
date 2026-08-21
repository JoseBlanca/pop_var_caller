//! Step 8 — how likely each genotype is **before any read is looked at**.
//!
//! At one locus, for one sample, this module produces one log-probability per candidate
//! genotype. The calling loop multiplies that by what the reads say and normalises; the
//! result is the posterior the caller emits (`doc/devel/ng/spec/calling_priors.md` §1).
//!
//! ## The belief has two sources, and both are measured before calling starts
//!
//! **How variable the population is.** In a nearly-invariant population almost every
//! sample carries two reference copies, and a single read showing something else is more
//! likely to be an error than a variant. In a diverse one it is not.
//!
//! **How homozygous this individual runs.** A selfing tomato accession is homozygous
//! nearly everywhere; an outbred human is not. That is the inbreeding coefficient, a
//! property of the sample rather than of the locus.
//!
//! Both arrive frozen from the parameter pre-pass — or, where the pre-pass could fit
//! nothing at all, the diversity is a species-range guess the run has to report as such
//! ([`ExpectedHeterozygosity::SPECIES_FALLBACK`](crate::ng::types::ExpectedHeterozygosity::SPECIES_FALLBACK)).
//! **This module fits nothing.** The one
//! thing it learns — how common each allele is *at this locus* — it learns from the other
//! samples in the cohort while the loop runs, and the sample's own contribution is
//! subtracted out so its reads cannot arrive twice (spec §6).
//!
//! ## Two functions with one contract between them
//!
//! 1. **Build a concentration** — per sample, per pass: the run's starting point plus what
//!    the *other* samples showed here. One addition per allele.
//! 2. **Turn a concentration into a log-prior row** — one number per candidate genotype.
//!    Only the differences between them matter: the loop adds what the reads say and
//!    rescales the row to sum to one, so a constant shared by every entry cancels. This is
//!    the costly half — one `lgamma` per allele a genotype carries a copy of, plus one
//!    `logsumexp` per homozygous genotype (arch §3.2) — and the loop runs it once per
//!    sample per pass.
//!
//! A **concentration** is one positive number per allele, read as *chromosomes the prior
//! behaves as though it had already seen*. Their ratio is the frequency it expects, their
//! sum is how much conviction that is (spec §1). Reading them as chromosome counts is what
//! makes the cohort term obvious: observed allele copies are added straight onto them,
//! because they are the same unit.
//!
//! ## What lands where
//!
//! - [`dirichlet_multinomial`] — the log-prior row itself. The Dirichlet-multinomial is
//!   what you get by drawing the locus's allele frequencies from a Dirichlet and then
//!   drawing this sample's few allele copies from those frequencies. Because it *averages
//!   over* the unknown frequencies rather than fixing them at one estimate —
//!   marginalizing — it is the default here (plan steps B1–B2).
//! - [`seed_spectrum`] — the SNP/indel starting point, read off the pre-pass's fitted
//!   frequency spectrum (plan step D).
//! - [`seed_ssr`] — the STR starting point: mass falling off geometrically from the
//!   cohort's modal repeat count, totalling what the cohort's measured repeat diversity
//!   implies (plan step E). **STR** in prose, `ssr` in module paths, as everywhere in ng.
//! - [`plug_in`] — the comparator: Hardy–Weinberg at a single estimated frequency, kept
//!   only so the change the marginalized prior makes stays measurable (plan step F).
//!
//! Build order and step contracts: `doc/devel/ng/impl_plan/calling_prior.md`. The design
//! and every why: `doc/devel/ng/spec/calling_priors.md`; the types and this seam:
//! `doc/devel/ng/arch/calling_priors.md`.

pub mod dirichlet_multinomial;
pub mod plug_in;
pub mod seed_spectrum;
pub mod seed_ssr;
