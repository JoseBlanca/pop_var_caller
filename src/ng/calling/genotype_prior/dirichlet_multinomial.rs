//! The log-prior row: a handful of allele copies drawn from a population whose composition
//! is not known. That is the Dirichlet-multinomial — the locus's allele frequencies drawn
//! from a Dirichlet and averaged out, with this sample's copies drawn from what is left.
//!
//! **Empty until plan step B1** (`doc/devel/ng/impl_plan/calling_prior.md`), which ports
//! production's Dirichlet-multinomial primitive (`src/genetics.rs`) into a form that fills
//! a caller's slice instead of returning a fresh `Vec`, and step B2, which wraps it in the
//! two-branch inbreeding mixture as `MarginalizedDirichletPrior` — the shipping default of
//! this step's seam (`doc/devel/ng/spec/calling_priors.md` §3).
