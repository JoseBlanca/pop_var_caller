//! The SNP/indel starting point: the run's two concentration numbers — the chromosomes'
//! worth of prior belief attached to the reference allele and to the alternatives — read
//! off the pre-pass's fitted frequency spectrum.
//!
//! **Empty until plan step D** (`doc/devel/ng/impl_plan/calling_prior.md`). At an ordinary
//! site most alternative alleles are rare, so the one chromosome the reference records is
//! almost always the common one and the reference's number is the larger — but only just,
//! and *how much* larger is fitted rather than fixed
//! (`doc/devel/ng/spec/calling_priors.md` §4.1).
