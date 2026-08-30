//! **A called locus becomes a record** — the mapper the output stream runs, and the seam where
//! this module meets the rest of the run.
//!
//! # This module is written against an interface that does not exist yet
//!
//! The run will be a stream: variants come off the caller, pass through filters and through
//! mappers that attach genotypes and annotations, and are written. This is one of those
//! mappers — the last one, turning what the caller produced into what the file states.
//!
//! **The part being assumed is small and named.** [`LocusInference`] already carries the calls,
//! the fitted copies, the allele table and whether the loop converged, and this reads all four
//! from the real type. What no type carries today is per-sample evidence that has outlived the
//! locus: how many reads matched each allele, how many the written alleles do not explain, the
//! cohort-pooled mapping qualities, the reference base beside the span, the tract's motif, and —
//! the one that is a genuine open question — **whether a sample's reads said anything at all**.
//! Those arrive in [`LocusEvidenceForOutput`], which is this module's guess at the shape the
//! stream will hand over.
//!
//! **What is settled and what is provisional**, so the next person can tell them apart:
//!
//! - *Settled* — everything the record does with these inputs. The no-call rule (spec §7.1), the
//!   `AF` denominator, the padding rule, the field set. Those are design decisions with
//!   documents behind them and they will not move when the interface does.
//! - *Provisional* — the shape of [`LocusEvidenceForOutput`] and [`SampleEvidenceForOutput`],
//!   and in particular where [`SampleEvidenceForOutput::reads_were_uninformative`] is computed.
//!   §"The one bit that needs the loop" says what it costs to supply and why it cannot be
//!   recovered here.
//!
//! When the stream lands, expect to rewrite the *construction* of these two structs and nothing
//! below them.

use super::{
    FilterVerdict, MapqPool, PaddingBase, SampleCall, SampleColumn, SampleReadCounts,
    TractAnnotation, VcfRecord,
};
use crate::ng::calling::quality::artifact_correction::ArtifactPenalties;
use crate::ng::calling::{LocusInference, SampleGenotypeCall};
use crate::ng::types::Phred;

/// **What one sample showed at a locus, once the locus itself is gone.**
///
/// Summed in the worker while the merge's `SampleSupport` is still in hand — the counts cannot
/// be recovered afterwards, which is why they travel rather than being re-derived.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct SampleEvidenceForOutput {
    /// Reads whose complete observation matched each allele of the locus's table, in allele
    /// order, reference first. Becomes `AD`.
    pub allele_reads: Vec<u32>,
    /// Reads this sample observed that no *written* allele explains: a dropped candidate's
    /// reads, and partial observations. Becomes `DP − ΣAD`.
    pub reads_no_written_allele_explains: u32,
    /// **Whether this sample's reads said nothing about which genotype it has** — its genotype
    /// likelihoods were flat.
    ///
    /// This is what turns a called sample into a `./.` (spec §7.1), and it is the one input
    /// here that **cannot be computed from anything downstream**. See the module's own note: the
    /// likelihoods live in per-sample scratch the loop overwrites, so whoever scores the sample
    /// has to answer this while it still knows, exactly as the genotype quality is taken during
    /// the loop's final pass rather than after it.
    ///
    /// **It must be the likelihood and not the posterior.** A sample with no reads is scored by
    /// the loop and comes back with a genotype, because the prior decides it alone — and at a
    /// locus where the fitted frequency is low that posterior is sharply peaked, so no threshold
    /// on it would catch this sample. The likelihood is what the reads said; the posterior is
    /// what the reads said plus what the cohort assumed.
    pub reads_were_uninformative: bool,
}

/// **Everything a record needs that the called locus does not already carry.**
///
/// Parallel to the locus's allele table where it says so, and checked against it by
/// [`assemble_record`] before anything is built.
#[derive(Clone, PartialEq, Debug)]
pub struct LocusEvidenceForOutput {
    /// One entry per sample of the run, in the run's sample order — the same order
    /// [`LocusInference::per_sample`] is in.
    pub samples: Vec<SampleEvidenceForOutput>,
    /// Cohort-pooled mapping qualities, one per allele, reference first.
    pub allele_mapq: Vec<MapqPool>,
    /// The reference base beside the span, resolved where the reference was still open —
    /// `Some` exactly when some allele of the locus is empty (spec §5).
    pub padding_base: Option<PaddingBase>,
    /// The site quality **after** the artifact correction, and the two penalties it subtracted.
    ///
    /// **Not read off the locus**, whose own quality field is the uncorrected baseline and which
    /// nothing between the worker and the correction stage may read
    /// (`doc/devel/ng/spec/calling_quality.md` §3.5).
    pub corrected_site_quality: Phred,
    /// What the two artifact tests took off, or `None` where they did not run.
    pub artifact_penalties: Option<ArtifactPenalties>,
    /// The tract's motif, at a repeat tract — `None` at a SNP or indel.
    ///
    /// **Not on the called locus today.** `RepeatTractProvenance` records which rung the tract's
    /// prior came from and how many scoring cells were defaulted, not the motif itself; the
    /// motif reaches the merge in `candidate_alleles_ssr.md`'s Milestone A, which is unbuilt.
    pub repeat_tract: Option<TractAnnotation>,
    /// The filter this locus is written on.
    ///
    /// **The emission steps' decision, not this module's** — it is a parameter for that reason.
    /// A locus whose loop did not converge is the one verdict derivable from the locus itself,
    /// and [`assemble_record`] checks the two agree rather than choosing between them.
    pub filter: FilterVerdict,
}

/// **Turn a called locus into the record the file will state.**
///
/// The mapper of the output stream. Every rule it applies is settled design; only the shape of
/// what it is handed is provisional (see the module note).
///
/// # Panics
///
/// When the evidence does not describe this locus: a sample list of a different length than the
/// locus's calls, a per-allele vector of the wrong width, or a filter that says the loop
/// converged when it did not. All three are wiring defects between two things the same worker
/// built, which is [`VcfRecord::new`]'s reasoning and the same choice.
#[must_use]
pub fn assemble_record(locus: &LocusInference, evidence: LocusEvidenceForOutput) -> VcfRecord {
    assert_eq!(
        evidence.samples.len(),
        locus.per_sample.len(),
        "the evidence carries {} samples and the locus was called on {}: both are the run's \
         sample order, so two lengths mean they were gathered over different cohorts",
        evidence.samples.len(),
        locus.per_sample.len()
    );
    let converged_filter = evidence.filter != FilterVerdict::EmDidNotConverge;
    assert_eq!(
        converged_filter,
        locus.converged,
        "the locus {} and its filter says {}: the loop's own answer is what that verdict \
         states, so the two disagreeing means one of them was set by hand",
        if locus.converged {
            "converged"
        } else {
            "did not converge"
        },
        evidence.filter.as_str()
    );

    let alleles: Vec<Box<[u8]>> = locus
        .alleles()
        .iter()
        .map(|bases| bases.to_vec().into_boxed_slice())
        .collect();

    let sample_columns = locus
        .per_sample
        .iter()
        .zip(evidence.samples)
        .map(|(call, sample)| SampleColumn {
            call: written_call(call, &sample),
            read_counts: SampleReadCounts::new(
                sample.allele_reads,
                sample.reads_no_written_allele_explains,
            ),
        })
        .collect();

    VcfRecord::new(
        locus.region,
        alleles,
        locus.cohort_expected_copies().copies().to_vec(),
        sample_columns,
        evidence.allele_mapq,
        evidence.padding_base,
        evidence.corrected_site_quality,
        evidence.artifact_penalties,
        evidence.filter,
        evidence.repeat_tract,
    )
}

/// **What the file says this sample is** — the no-call rule of spec §7.1, in one place.
///
/// Two routes to `./.`, and they are different facts about the sample:
///
/// - **the caller declined to invent a genotype** — candidate selection cut an allele this
///   sample's own reads had earned, so the locus is called over a set that cannot hold what it
///   carries ([`SampleGenotypeCall::Missing`]);
/// - **the sample's reads said nothing** — its likelihoods were flat, so the genotype the loop
///   produced came from the prior alone. Spec §7.1: a sample is never force-called for lack of
///   evidence.
///
/// A sample that is neither keeps the genotype and quality the loop gave it.
fn written_call(call: &SampleGenotypeCall, evidence: &SampleEvidenceForOutput) -> SampleCall {
    match call {
        SampleGenotypeCall::Missing => SampleCall::NoCall,
        SampleGenotypeCall::Called {
            genotype,
            genotype_quality,
        } => {
            if evidence.reads_were_uninformative {
                SampleCall::NoCall
            } else {
                SampleCall::Called {
                    genotype: genotype.clone(),
                    genotype_quality: *genotype_quality,
                }
            }
        }
    }
}

#[cfg(test)]
mod tests;
