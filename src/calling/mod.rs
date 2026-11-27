//! Turning a locus's reads into a call: which alleles are in play, how probable each
//! sample's reads are under each candidate genotype, how likely each genotype is before
//! the reads are looked at, and what comes out. Steps 6 to 9 of the ng caller.
//!
//! **One folder rather than four, and a dependency decided it.** As four sibling
//! folders at the top of `src/`, the three calling architecture documents each had
//! to state by hand that the prior and the likelihood must never import from the
//! inference step — a rule the tree made necessary and could not enforce. Under one
//! folder the types all four steps share sit one level above them, in the files
//! listed below, and each step imports them from there
//! (`doc/devel/ng/arch/module_layout.md`, principle 1b). The rule still holds — no code
//! in [`genotype_prior`] or [`likelihood`] imports from [`inference`] — though nothing
//! but review keeps it so.
//!
//! **The shared types are not free of the sub-modules, and imports run both ways
//! between the two levels.** Only [`CandidateAlleles`] and [`ExpectedAlleleCopies`]
//! need nothing from below. [`LocusEvidence`] carries the likelihood's per-sample
//! evidence views; [`FrozenParameters`] the likelihood's calibration and contamination
//! and the prior's seed; [`LocusInference`] the site-quality step's artifact counts —
//! in each case a type the step that produces it defines. [`CallingScratch`] goes
//! further by design: it holds a buffer from every sub-module, the inference step's
//! included, because a worker owns one scratch and lends each step its part.
//!
//! Steps 10 and 11 are deliberately elsewhere: phasing and locus filtering read what
//! calling produces rather than taking part in producing it.
//!
//! ## What is here
//!
//! The shared types live one family to a file, and all are re-exported here, so each is
//! named `crate::calling::X`:
//!
//! - `candidate_alleles.rs` — [`CandidateAlleles`], the alleles a locus is called over,
//!   and [`ExpectedAlleleCopies`], the fractional allele counts the loop feeds back to
//!   itself.
//! - `locus_evidence.rs` — [`LocusEvidence`], one locus's reads per sample and the one
//!   place the SNP/indel and repeat-tract paths meet, with [`GenericLocusSample`] and
//!   [`SsrLocusShape`].
//! - `frozen_parameters.rs` — [`FrozenParameters`], everything the parameter pre-pass
//!   fitted, gathered so one borrow crosses the seam.
//! - `scratch.rs` — [`CallingScratch`], every buffer a locus's calling fills, allocated
//!   once per worker, and the borrow structs that hand each step its share of it.
//! - `locus_inference.rs` — [`LocusInference`] with [`SampleGenotypeCall`] and
//!   [`RepeatTractProvenance`], what a locus produces.
//! - `tests.rs`, test-only — the fixtures the tests of more than one family use, and
//!   that `inference`'s tests borrow.
//!
//! Of the four sub-modules — the candidate
//! step, the likelihood, the genotype prior and the inference loop — three are here:
//! [`allele_candidates`], step 6, which narrows the merge's allele table to the sequences
//! worth calling over (`doc/devel/ng/impl_plan/candidate_alleles.md`);
//! [`likelihood`], step 7, how probable this sample's reads are given each genotype
//! (`doc/devel/ng/impl_plan/calling_read_likelihoods.md`); and
//! [`genotype_prior`], step 8, how likely each genotype is before any read is looked at
//! (`doc/devel/ng/impl_plan/calling_prior.md`). **The fourth is [`inference`], step 9, the
//! loop that consumes all three** — it holds the seam every way of handling a cohort crosses
//! and the configuration of the three nested loops, two of which ship switched off
//! (`doc/devel/ng/impl_plan/calling_loop.md`).
//!
//! Beside this file rather than inside any one sub-module: [`genotype_table`], which
//! says which genotypes a locus's alleles make and holds the three flat tables the
//! genotype prior and the read likelihood both read them from. It sits here for the
//! same reason the types above do — two of the four sub-modules consume it, so it
//! belongs one level above them — and its three types are re-exported here, so that
//! all of calling's shared vocabulary is named at one depth. Beside it again, and
//! test-only: `genotype_table_parity`, which checks that table value for value against
//! production's `GenotypeShape`, frozen in `testdata/` at promotion step C5. It is a file of its
//! own rather than a block inside [`genotype_table`]'s own tests so that the table's tests stay
//! free of the comparison.
//!
//! The scalars other steps also name — `AlleleId`, `Phred`, `Genotype` — are not here
//! but in [`crate::types`], with the rest of ng's shared vocabulary.

#[cfg(test)]
mod genotype_table_parity;
#[cfg(test)]
mod loop_parity;
#[cfg(test)]
mod quality_parity;

pub mod allele_candidates;
pub mod evidence_shaping;
pub mod genotype_prior;
pub mod genotype_table;
pub mod inference;
pub mod likelihood;
pub mod parameters_file;
pub mod quality;
pub mod run_parameters;

mod candidate_alleles;
mod frozen_parameters;
mod locus_evidence;
mod locus_inference;
mod scratch;

pub use candidate_alleles::{CandidateAlleles, ExpectedAlleleCopies};
pub use frozen_parameters::FrozenParameters;
pub use genotype_table::{GenotypeIdx, GenotypeTable, GenotypeTableView};
pub use likelihood::generic::{
    ERROR_SPREAD_BASES, NO_ERROR_SPREAD, ReferenceBesideLocus, WrongReadCharge,
    WrongReadChargeTable, allele_is_compatible_with_partial, fill_wrong_read_charges,
    partial_row_fits_allele,
};
pub use likelihood::{
    ContaminationMixture, ContaminationView, GenericObservation, GenericRowScratch,
    GenericSampleEvidence, MIN_BASE_ERROR, MIN_CONTAMINANT_FREQUENCY, ReadGroupCalibration,
    ReadGroupParameters, SsrRowScratch, SsrSampleEvidence, fill_batch_allele_copies,
    fill_own_batch_contaminant_row, fill_shared_contaminant_allele_frequencies,
};
pub use locus_evidence::{GenericLocusSample, LocusEvidence, SsrLocusShape};
#[cfg(test)]
pub use locus_inference::placeholder_genotype_rows;
pub use locus_inference::{LocusInference, RepeatTractProvenance, SampleGenotypeCall};
// `EmissionCost` is named only by tests and by doc links in `inference`, which the unused-import
// lint does not count.
#[allow(unused_imports)]
pub(crate) use scratch::EmissionCost;
pub use scratch::{CallingScratch, UNWRITTEN_SCRATCH_VALUE};
pub(crate) use scratch::{CohortSummingBuffers, SampleScoringBuffers};

#[cfg(test)]
mod tests;
