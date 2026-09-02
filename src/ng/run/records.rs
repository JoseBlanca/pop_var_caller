//! **What a record needs that the called locus does not carry** — the run's half of the seam
//! [`ng::vcf::assemble`](crate::ng::vcf::assemble) describes from the other side.
//!
//! `assemble_record` turns a [`LocusInference`] plus a [`LocusEvidenceForOutput`] into a
//! record, and its module note says the second of those is *"this module's guess at the shape
//! the stream will hand over"*. This is the stream. Everything here is read while the merge's
//! [`CohortObservation`] is still in hand, because none of it can be recovered afterwards: the
//! observation is dropped as soon as its genotypes exist, which is what calling inside the
//! builder buys (`doc/devel/ng/arch/run_streaming.md` §3.4).
//!
//! # The three numberings meet here again, and one of them is new
//!
//! [`shape_generic_locus`](crate::ng::calling::evidence_shaping) already has to reconcile the
//! merge's covering samples, the run's sample order and the calling scratch's rows. The record
//! adds a fourth axis that is not a sample axis at all: **the merge's allele table and the
//! record's are different tables**, because candidate selection drops alleles between them. So
//! every per-allele count here is keyed by [`AlleleRemap`] and never by the merge's own index —
//! a read on an allele selection dropped belongs in `DP` and in no `AD` slot
//! (`doc/devel/ng/spec/vcf_output.md` §7), and using the merge's index would put it in the
//! wrong one instead.
//!
//! # The padding base, and why it is fetched here
//!
//! VCF cannot spell an empty allele, so a record with one is written by giving every allele a
//! flanking reference base (spec §5). **The locus does not carry that base**: the merge gathers
//! the bases of the span and no more. So it is fetched from the run's own reference — one
//! accessor for the whole run, minted beside the walkers' — at the moment the record is built,
//! which is the only moment that knows both the span and whether any allele is empty.
//!
//! **⚑ And on today's path it is never fetched, which is worth knowing before reading this as
//! live code.** The generic mint anchors its indels: an insertion's reference span is its anchor
//! base alone and a deletion's is the anchor plus the deleted run
//! (`ReadEvent::footprint_span`, `locus_generation/pileup/decompose.rs`), so a
//! deletion's alternative is one base and never nothing, and **no allele a generic locus is
//! called over is ever empty**. The empty allele spec §5 was written for is the repeat-tract
//! path's full-tract deletion, and that path is unbuilt. What this is, then, is the answer
//! ready for when a record needs it — and it cannot be left out: `VcfRecord::new` asserts a
//! padding base is carried **exactly** when some allele is empty, so a run that did not compute
//! one would panic at the first record that needed it rather than write a wrong one.
//!
//! **Never invented.** Production's repeat-tract writer puts the letter `N` at a span starting
//! at a contig's first base ([`vcf_out.rs:405-435`](../../../../src/ssr/cohort/vcf_out.rs));
//! spec §5 declines to port that, so a base that cannot be read is a run that stops and says
//! which position it could not read.

use crate::ng::calling::allele_candidates::{AlleleRemap, SelectionVerdict, UnmatchedSupport};
use crate::ng::calling::quality::artifact_correction::correct_site_quality;
use crate::ng::calling::{CandidateAlleles, LocusInference, SampleGenotypeCall};
use crate::ng::locus_generation::LocusKind;
use crate::ng::ref_seq::{EvictableRefSeq, RefSeq, RefSeqError};
use crate::ng::run::cohort_merge::build::CohortObservation;
use crate::ng::types::{GenomeRegion, Position};
use crate::ng::vcf::assemble::{LocusEvidenceForOutput, SampleEvidenceForOutput};
use crate::ng::vcf::{FilterVerdict, MapqPool, PaddingBase, TractAnnotation};

/// The first base of every contig, 1-based — the one position with no base to its left.
const FIRST_POSITION_OF_A_CONTIG: Position = Position(1);

/// **Whether this locus establishes a variant at all**, which is what decides if it reaches
/// the file (`doc/devel/ng/spec/vcf_output.md` §9).
///
/// A locus every sample was called homozygous-reference at, or no sample was called at,
/// established nothing: its absence from the file says *nothing here*, which is exactly what
/// production's generic path means by leaving it out. There is no gVCF and no reference block.
///
/// **The test is on the calls the file would write, not on the ones the loop made.** A sample
/// whose reads said nothing is written `./.` (spec §7.1) and carries no allele, so it cannot be
/// what keeps a locus in the file — which is the same rule
/// [`assemble_record`](crate::ng::vcf::assemble::assemble_record) applies one step later.
#[must_use]
pub fn a_written_genotype_carries_an_alternative(locus: &LocusInference) -> bool {
    locus.per_sample.iter().any(|call| match call {
        SampleGenotypeCall::Missing => false,
        SampleGenotypeCall::Called {
            genotype,
            reads_were_uninformative,
            ..
        } => {
            !*reads_were_uninformative
                && genotype
                    .alleles()
                    .iter()
                    .any(|allele| !allele.is_reference())
        }
    })
}

/// **The reference base a record with an empty allele is padded with**, or `None` where every
/// allele spells bases.
///
/// Ordinarily the base immediately to the **left** of the span, with `POS` moving one base left
/// with it. At a span starting at the contig's first base there is nothing to the left, and the
/// base immediately to the **right** is appended instead with `POS` unmoved — the VCF 4.4 rule
/// for an event at position 1 (spec §5).
///
/// `scratch` is the fetch's destination and is cleared by the fetch; it is a parameter so that
/// a run reuses one buffer over a genome rather than allocating a one-byte `Vec` per indel.
///
/// # Errors
///
/// Whatever the reference fetch refuses. The ordinary reachable cause is the FASTA becoming
/// unreadable part-way through a run; the position being out of bounds needs a span that
/// starts at a contig's first base *and* ends at its last, which no aligned read can produce —
/// a read carrying a deletion must match reference on at least one side of it.
pub fn padding_base_beside<R>(
    reference: &R,
    region: GenomeRegion,
    alleles: &CandidateAlleles,
    scratch: &mut Vec<u8>,
) -> Result<Option<PaddingBase>, RefSeqError>
where
    R: RefSeq + EvictableRefSeq,
{
    if !alleles.iter().any(<[u8]>::is_empty) {
        return Ok(None);
    }
    let at_contig_start = region.start == FIRST_POSITION_OF_A_CONTIG;
    let position = if at_contig_start {
        region.end.get() + 1
    } else {
        region.start.get() - 1
    };
    reference.fetch_into(region.contig, position, 1, scratch)?;
    let base = scratch[0];
    // **Release what the run has walked past**, so the accessor's window does not grow into a
    // resident contig over a genome. Correctness does not depend on it: an evicted position is
    // simply read again, and the loci this is driven over arrive in genome order, so nothing
    // asks for a base behind the one just fetched.
    reference.evict_before(position);
    Ok(Some(if at_contig_start {
        PaddingBase::Right(base)
    } else {
        PaddingBase::Left(base)
    }))
}

/// **Everything a record needs that the called locus does not already carry**, gathered while
/// the cohort observation is still in hand.
///
/// The four kinds of thing it collects, and where each comes from:
///
/// - **what each sample's reads showed** — `AD` per written allele, and how many of the
///   sample's reads no written allele explains, which is `DP − ΣAD` (spec §7). Read off the
///   merge's per-sample support rows through `remap`, and off what candidate selection set
///   aside;
/// - **the cohort's mapping qualities per written allele** — the same reads `AD` counts, which
///   is why they are summed from the same rows;
/// - **the site quality after the artifact correction**, and the two penalties it charged.
///   `LocusInference`'s own quality field is the **uncorrected** baseline and nothing between
///   the worker and this correction may read it as a site quality
///   (`doc/devel/ng/spec/calling_quality.md` §3.5);
/// - **the padding base**, already resolved by [`padding_base_beside`].
///
/// # Panics
///
/// On a covering sample naming a run sample the inference was not called over — the merge and
/// the calling loop are handed the same run sample count, so the two disagreeing is a wiring
/// defect in the driver that built both, which is
/// [`assemble_record`](crate::ng::vcf::assemble::assemble_record)'s reasoning for its own
/// checks and the same choice.
#[must_use]
pub fn evidence_for_output(
    locus: &LocusInference,
    observation: &CohortObservation,
    remap: &AlleleRemap,
    unmatched: &[UnmatchedSupport],
    selection: SelectionVerdict,
    padding_base: Option<PaddingBase>,
) -> LocusEvidenceForOutput {
    let written_alleles = locus.alleles().len();
    let run_sample_count = locus.per_sample.len();
    let mut samples: Vec<SampleEvidenceForOutput> = (0..run_sample_count)
        .map(|_| SampleEvidenceForOutput {
            allele_reads: vec![0; written_alleles],
            reads_no_written_allele_explains: 0,
        })
        .collect();
    let mut allele_mapq = vec![
        MapqPool {
            reads: 0,
            mapq_sum: 0
        };
        written_alleles
    ];

    for (covering, support) in observation.per_sample.iter().enumerate() {
        assert!(
            support.sample < run_sample_count,
            "the merge's entry {covering} at {} names run sample {} and the locus was called \
             over {run_sample_count}: both are the run's sample order, so an index past its end \
             means the two were built over different cohorts",
            observation.region,
            support.sample
        );
        let sample = &mut samples[support.sample];
        for row in &support.supported {
            // **The merge's allele index is not the record's.** A row whose allele candidate
            // selection dropped has no slot to be counted in, and its reads reach `DP` through
            // the leftover below rather than through any `AD`.
            let Some(written) = remap.candidate_for(row.allele) else {
                continue;
            };
            let written = usize::from(written.get());
            sample.allele_reads[written] =
                sample.allele_reads[written].saturating_add(row.support.num_reads);
            allele_mapq[written].reads += u64::from(row.support.num_reads);
            allele_mapq[written].mapq_sum += u64::from(row.support.mapq_sum);
        }
        // **Three kinds of read no written allele explains, and they are disjoint by
        // construction.** The leftover counts this sample's reads on sequences selection
        // dropped — a scan of the same support rows, asked by the step that dropped them. A
        // partial observation never reached the allele table at all: it says the sample carries
        // *at least* this much and is scored on its own axis. And a read removed as evidence —
        // named at some of this sample's records inside the locus and not at all of them —
        // reaches no `supported` row either, which is why the merge calls it lost depth and
        // counts it rather than leaving it to be inferred from an absence.
        //
        // **⚑ The third is a reading of spec §7, taken 2026-09-01, and the owner may overturn
        // it.** `DP` is *"every read observation the sample had at the locus, whether or not a
        // written allele explains them"*, and these reads were observed there: leaving them out
        // makes `DP` understate the depth at exactly the loci that span several of a sample's
        // records. `SampleEvidenceForOutput`'s own doc names two of the three, written before
        // this seam existed and calling its own shape provisional; `UnmatchedSupport`'s doc
        // excludes them from the *leftover* for a different reason — they carry no quality sum
        // and a read likelihood needs one — which settles the likelihood and not the file.
        // Disjointness is what makes adding them safe: the leftover is computed from
        // `supported`, and nothing these reads showed reaches it.
        let dropped = unmatched
            .get(covering)
            .map_or(0, |leftover| leftover.num_reads);
        let partial = support.partials.iter().fold(0_u32, |total, partial| {
            total.saturating_add(partial.num_reads)
        });
        let partial = partial.saturating_add(support.reads_removed_as_evidence);
        // **Added into the row rather than assigned over it**, so this count accumulates the
        // way `allele_reads` above does. The merge holds one entry per covering sample, so
        // nothing reaches a row twice; an assignment would nonetheless keep only the last of
        // two, which is a different answer from the one the alleles beside it would give.
        sample.reads_no_written_allele_explains = sample
            .reads_no_written_allele_explains
            .saturating_add(dropped.saturating_add(partial));
    }

    let (corrected_site_quality, artifact_penalties) = match locus.artifact_test_counts() {
        Some(counts) => {
            let (quality, penalties) =
                correct_site_quality(locus.uncorrected_site_quality(), &counts);
            (quality, Some(penalties))
        }
        // **A locus that gave the two tests nothing to weigh keeps its baseline**, and the
        // record says the tests did not run rather than writing two zeroed penalties, which a
        // reader could not tell from two tests that charged nothing.
        None => (locus.uncorrected_site_quality(), None),
    };

    LocusEvidenceForOutput {
        samples,
        allele_mapq,
        padding_base,
        corrected_site_quality,
        artifact_penalties,
        // **Read off the called locus's own candidate table**, which is where selection
        // stamped the kind: a repeat tract's table is `LocusKind::Ssr` and carries the motif,
        // and everything the record says about the repeat — the `STR` flag, `RU`, `PERIOD`
        // and each called allele's `REPCN` — is written from that one motif.
        repeat_tract: match locus.alleles().kind() {
            LocusKind::Ssr(detail) => Some(TractAnnotation::new(detail.motif)),
            LocusKind::Generic | LocusKind::SsrBundle => None,
        },
        filter: filter_for(locus, selection),
    }
}

/// **The one verdict this record is written on**, from the loop's answer and selection's.
///
/// # A loop that did not settle outranks everything selection could say
///
/// Not a preference: [`assemble_record`](crate::ng::vcf::assemble::assemble_record) asserts that
/// the filter is `EMNoConv` exactly when the loop failed to converge, so the two cannot be
/// reported together and the loop's answer is the one the file has to carry. A tract that was
/// truncated *and* did not converge says `EMNoConv`, and the truncation is visible in the record
/// anyway — the alternatives it kept are the alternatives it kept.
///
/// # The tract verdicts are a tract's, and one of them ng never mints
///
/// `notPeriodic` and `tooManyAlleles` come from repeat-tract selection
/// (`doc/devel/ng/spec/candidate_alleles_ssr.md` §6, §7). **`lowDepth` is declared and never
/// written**: production refuses a tract whose *cohort-summed* depth is under ten, and ng does
/// not port that gate — depth is asked once, upstream, per sample, by the merge's keep rule,
/// and there is no depth verdict on this path (§6). The vocabulary stays in the header because
/// a file written by an older caller can carry it.
///
/// # A truncated SNP/indel locus is not filtered, and that is unchanged
///
/// The ordinary path's cap cuts the lowest-ranked alternatives and calls the locus over the
/// rest; `tooManyAlleles` is spec §8's *tract* filter and stays one. What changed here is only
/// that a tract's truncation now reaches the column.
fn filter_for(locus: &LocusInference, selection: SelectionVerdict) -> FilterVerdict {
    if !locus.converged {
        return FilterVerdict::EmDidNotConverge;
    }
    match locus.alleles().kind() {
        LocusKind::Generic | LocusKind::SsrBundle => FilterVerdict::Pass,
        LocusKind::Ssr(_) => match selection {
            SelectionVerdict::NotPeriodic => FilterVerdict::NotPeriodic,
            SelectionVerdict::Truncated { .. } => FilterVerdict::TooManyAlleles,
            _ => FilterVerdict::Pass,
        },
    }
}

#[cfg(test)]
mod tests;
