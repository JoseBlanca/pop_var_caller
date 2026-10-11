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
//! - *Provisional* — the shape of [`LocusEvidenceForOutput`] and [`SampleEvidenceForOutput`].
//!
//! **⛦ The one input this module could not be given now comes from the loop (2026-09-01).**
//! Whether a sample's reads said anything is `SampleGenotypeCall::Called`'s
//! `reads_were_uninformative`, minted where the sample is scored — the likelihoods live in
//! per-sample scratch the next locus overwrites, so nowhere else can answer it. It used to be a
//! field of [`SampleEvidenceForOutput`] and this module's one open question.
//!
//! When the stream lands, expect to rewrite the *construction* of these two structs and nothing
//! below them.

use super::TractWidening;
use super::{
    FilterVerdict, MapqPool, PaddingBase, SampleCall, SampleColumn, SampleReadCounts,
    TractAnnotation, VcfRecord,
};
use crate::calling::quality::artifact_correction::ArtifactPenalties;
use crate::calling::{LocusInference, SampleGenotypeCall};
use crate::locus_generation::{FlankDeletion, LocusKind};
use crate::region_typing::run_parse::parse_into_runs;
use crate::types::{AlleleId, Genotype, Phred};
use crate::types::{GenomeRegion, Position};
use crate::window_coverage::WindowCoverage;

/// **What one sample showed at a locus, once the locus itself is gone.**
///
/// Summed in the worker while the merge's `SampleSupport` is still in hand — the counts cannot
/// be recovered afterwards, which is why they travel rather than being re-derived.
///
/// **Comparing two of these compares the window by bit pattern**, which is
/// [`WindowCoverage`]'s own `PartialEq` and not the float `==` a reader might expect: it is what
/// makes two absent windows — two pairs of `NaN`s — equal, and what would report a correct run as
/// disagreeing at every window the position floor silenced if it were `==` instead.
///
/// **`Eq` is gone from the derive** only because a derive cannot reach past a field whose type
/// withholds it, and `WindowCoverage` withholds it deliberately: bitwise equality separates `+0.0`
/// from `-0.0`, which is what comparing two runs wants and what a lookup key must not have. It was
/// dead weight here in any case — nothing in the tree compares two of these, and
/// [`LocusEvidenceForOutput`], the only type that holds them, never derived `Eq` either.
#[derive(Clone, PartialEq, Debug)]
pub struct SampleEvidenceForOutput {
    /// Reads whose complete observation matched each allele of the locus's table, in allele
    /// order, reference first. Becomes `AD`.
    pub allele_reads: Vec<u32>,
    /// Reads this sample observed that no *written* allele explains: a dropped candidate's
    /// reads, and partial observations. Becomes `DP − ΣAD`.
    pub reads_no_written_allele_explains: u32,
    /// This sample's read depth and GC over the 500-base window centred on the locus's first
    /// base (`doc/devel/ng/spec/window_coverage.md` §3.5), **or an absent pair where it has no
    /// usable window there**.
    ///
    /// **Dense over the run's samples, where the merge's is over the covering ones**: a sample
    /// that covered nothing at this locus has a row here and no entry there, and both say the
    /// same thing about its window — that there is none. It reaches no VCF field today; the
    /// hidden-duplication filter is what reads it.
    pub window_coverage: WindowCoverage,
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
    /// `Some` exactly when some allele of the locus is empty (spec §5). [`assemble_record`] drops
    /// it again when the only empty allele is an alternative no sample calls.
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

    let sample_columns: Vec<SampleColumn> = locus
        .per_sample
        .iter()
        .zip(evidence.samples)
        .map(|(call, sample)| SampleColumn {
            call: written_call(call),
            read_counts: SampleReadCounts::new(
                sample.allele_reads,
                sample.reads_no_written_allele_explains,
            ),
        })
        .collect();

    let mut expected_copies = locus.cohort_expected_copies().copies().to_vec();
    let mut allele_mapq = evidence.allele_mapq;
    let mut alleles = alleles;
    let mut sample_columns = sample_columns;
    let mut flank_deletions: Vec<FlankDeletion> = (0..alleles.len())
        .map(|index| locus.alleles().flank_deletion(index))
        .collect();
    drop_alternatives_no_sample_calls(
        &mut alleles,
        &mut expected_copies,
        &mut allele_mapq,
        &mut sample_columns,
        &mut flank_deletions,
    );
    let mut region = locus.region;
    let mut repeat_tract = evidence.repeat_tract;
    let mut padding_base = evidence.padding_base;
    // **An STR bundle's run counts are read off the alleles the record keeps**, before any
    // widening adds flank bases to them: each is parsed against the bundle's reference runs, the
    // reference allele being the reference hull (`str_bundle_loci.md` §7.1, §7.5).
    if let LocusKind::SsrBundle(detail) = locus.alleles().kind() {
        let counts = alleles
            .iter()
            .map(|allele| {
                parse_into_runs(allele, &alleles[0], &detail.runs)
                    .copies()
                    .collect()
            })
            .collect();
        repeat_tract = repeat_tract.map(|tract| tract.with_run_counts_of_alleles(counts));
    }
    let flanks = match locus.alleles().kind() {
        LocusKind::Ssr(detail) => Some((&*detail.left_flank, &*detail.right_flank)),
        LocusKind::SsrBundle(detail) => Some((&*detail.left_flank, &*detail.right_flank)),
        LocusKind::Generic => None,
    };
    if let Some(flanks) = flanks
        && let Some(widened) = WidenedTract::over(&alleles, &flank_deletions, flanks)
    {
        region = GenomeRegion {
            contig: region.contig,
            start: Position(region.start.get() - widened.left as u64),
            end: Position(region.end.get() + widened.right as u64),
        };
        repeat_tract = repeat_tract.map(|tract| {
            tract.widened(TractWidening {
                left: u8::try_from(widened.left).expect("a flank holds fewer than 256 bases"),
                right: u8::try_from(widened.right).expect("a flank holds fewer than 256 bases"),
                flank_bases_per_allele: widened.flank_bases_per_allele,
                reach_per_allele: widened.reach_per_allele,
            })
        });
        padding_base = widened.padding_base;
        alleles = widened.alleles;
    }
    // **The padding base was resolved over the whole candidate table, and the drop above can
    // take away the only empty allele** — a deletion that survived selection but that no sample
    // is called carrying. A record with no empty allele is written without padding (spec §5),
    // so the base goes with the allele that needed it. The drop removes alleles and never
    // shortens one, so it cannot make a padding base necessary where there was none.
    let padding_base = padding_base.filter(|_| alleles.iter().any(|allele| allele.is_empty()));

    VcfRecord::new(
        region,
        alleles,
        expected_copies,
        sample_columns,
        allele_mapq,
        padding_base,
        evidence.corrected_site_quality,
        evidence.artifact_penalties,
        evidence.filter,
        repeat_tract,
    )
}

/// **An STR record's alleles respelled over a reference widened into its flanks**
/// (`str_edge_deletion_alleles.md` §5.8) — `None` when no allele reaches into a flank, and the
/// record is written over the tract as before.
struct WidenedTract {
    /// Flank bases added before the tract, and after it.
    left: usize,
    right: usize,
    /// Every allele over the widened reference, the reference first.
    alleles: Vec<Box<[u8]>>,
    /// How many flank bases each allele was given around its spelling, for `REPCN`.
    flank_bases_per_allele: Vec<u8>,
    /// Each allele's reach into the left and the right flank, as the generic path places it.
    reach_per_allele: Vec<(u8, u8)>,
    /// The base before the widened record, for when an allele comes out empty.
    padding_base: Option<PaddingBase>,
}

impl WidenedTract {
    /// Widen over each allele's reach into the flanks as the generic path would place it.
    ///
    /// **Right**: an allele's flank deletion's right count — its smallest form already places
    /// a deletion as far left as it goes on that side. **Left**: how many left-flank bases,
    /// counted back from the tract, lie at or after the first base the generic path shows as
    /// changed. The allele's sequence there is the left flank less its deleted bases, then its
    /// spelling, against the reference's left flank and tract. The generic path left-aligns
    /// every read, so a deletion starts as far left as its sequence allows: where the two
    /// sequences' shared ending runs out, less the deletion's length — or, for an allele that
    /// is not shorter, where they first differ. For a plain allele the reach is zero unless
    /// its deletion slides into the flank, as `AG` from `…CG|AGAG…` does.
    fn over(
        alleles: &[Box<[u8]>],
        flank_deletions: &[FlankDeletion],
        (left_flank, right_flank): (&[u8], &[u8]),
    ) -> Option<Self> {
        let flank = left_flank.len();
        let reference = &alleles[0];
        let reach_left = |bases: &[u8], deleted: FlankDeletion| -> usize {
            let kept = flank.saturating_sub(usize::from(deleted.left));
            let ours: Vec<u8> = left_flank[..kept].iter().chain(bases).copied().collect();
            // The reference over the same stretch: the allele also covers the right-flank bases
            // its extent takes there.
            let right_extent = usize::from(deleted.right).min(right_flank.len());
            let theirs: Vec<u8> = left_flank
                .iter()
                .chain(reference.iter())
                .chain(&right_flank[..right_extent])
                .copied()
                .collect();
            let first_difference = ours.iter().zip(&theirs).take_while(|(a, b)| a == b).count();
            let start = if ours.len() < theirs.len() {
                let shared_ending = ours
                    .iter()
                    .rev()
                    .zip(theirs.iter().rev())
                    .take_while(|(a, b)| a == b)
                    .count();
                first_difference.min(theirs.len() - shared_ending - (theirs.len() - ours.len()))
            } else {
                first_difference
            };
            flank.saturating_sub(start)
        };
        let reach_per_allele: Vec<(usize, usize)> = alleles
            .iter()
            .zip(flank_deletions)
            .enumerate()
            .map(|(index, (bases, deleted))| {
                if index == 0 {
                    (0, 0)
                } else {
                    (reach_left(bases, *deleted), usize::from(deleted.right))
                }
            })
            .collect();
        let left = reach_per_allele
            .iter()
            .map(|reach| reach.0)
            .max()
            .unwrap_or(0);
        let right = reach_per_allele
            .iter()
            .map(|reach| reach.1)
            .max()
            .unwrap_or(0);
        if left == 0 && right == 0 {
            return None;
        }
        // A reach past the flank the locus holds cannot be written over it; read measurement
        // records no such deletion, so this is a plain allele sliding the whole flank's length.
        if left > flank || right > right_flank.len() {
            return None;
        }
        let mut widened = Vec::with_capacity(alleles.len());
        let mut flank_bases_per_allele = Vec::with_capacity(alleles.len());
        for (bases, deleted) in alleles.iter().zip(flank_deletions) {
            let (cut_left, cut_right) = (usize::from(deleted.left), usize::from(deleted.right));
            let before = &left_flank[flank - left..flank - cut_left.min(left)];
            let after = &right_flank[cut_right.min(right)..right];
            let mut allele = Vec::with_capacity(before.len() + bases.len() + after.len());
            allele.extend_from_slice(before);
            allele.extend_from_slice(bases);
            allele.extend_from_slice(after);
            flank_bases_per_allele.push(
                u8::try_from(before.len() + after.len())
                    .expect("two flanks of a tract hold fewer than 256 bases"),
            );
            widened.push(allele.into_boxed_slice());
        }
        // An allele that comes out empty is padded with the base before the widened record,
        // or the one after it where the flank is spent.
        let padding_base = if left < flank {
            Some(PaddingBase::Left(left_flank[flank - left - 1]))
        } else if right < right_flank.len() {
            Some(PaddingBase::Right(right_flank[right]))
        } else {
            None
        };
        if padding_base.is_none() && widened.iter().any(|allele| allele.is_empty()) {
            return None;
        }
        Some(Self {
            left,
            right,
            alleles: widened,
            flank_bases_per_allele,
            reach_per_allele: reach_per_allele
                .into_iter()
                .map(|(l, r)| {
                    (
                        u8::try_from(l).expect("a flank holds fewer than 256 bases"),
                        u8::try_from(r).expect("a flank holds fewer than 256 bases"),
                    )
                })
                .collect(),
            padding_base,
        })
    }
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
/// **Drop every alternative allele no sample's genotype names**, and renumber what is left.
///
/// # What it removes, and why it is not "alleles with no evidence"
///
/// A record's alternatives are the candidates the locus was *called over*, not the ones the
/// cohort turned out to carry. A candidate can survive selection on a couple of reads, take no
/// called copy in any sample, and still reach the file — measured on
/// `benchmarks/giab/per_sample` at 30×, 9 records of 2,290 carried one, and each such allele had
/// two reads behind it. **So the test cannot be "has no reads", and it cannot be applied where
/// the observations are made**: at that point nothing knows yet what will be called. It is
/// `AC == 0` — no genotype in the file names this allele — and it can only be asked once every
/// sample has been genotyped, which is here.
///
/// It is what freebayes does (its `MEANALT` says how many alternatives it saw at a site where it
/// reports one) and what `bcftools view --trim-alt-alleles` does to a file after the fact.
///
/// # The reads do not vanish; they move one column
///
/// A dropped allele's `AD` entry is added to that sample's `DP − ΣAD`, which ng already writes
/// and already documents as *reads no written allele explains*
/// ([`SampleReadCounts`](crate::vcf::SampleReadCounts)) — the column that already holds a
/// dropped candidate's reads and the partial observations. So the depth is conserved and lands
/// in the field whose meaning already covers it.
///
/// # When no sample names anything, nothing is dropped
///
/// **A refused locus writes every sample as a no-call** (spec §8), so no genotype names any
/// allele and this rule, applied literally, would empty the ALT column of a record whose whole
/// purpose is to say the caller could not decide. The same holds for a cohort genotyped
/// homozygous-reference throughout. In both cases the record is left exactly as it was: the
/// alternatives it carries are the only statement it has to make.
///
/// # What is left alone
///
/// `QUAL`, the site's artifact penalties and the filter verdict were computed over the whole
/// candidate set and are **not** recomputed. That is the same choice `bcftools` makes, and the
/// honest reading of the field: the quality is the evidence the site was judged on, not a
/// property of the alleles that survived the judging.
pub(crate) fn drop_alternatives_no_sample_calls(
    alleles: &mut Vec<Box<[u8]>>,
    expected_copies: &mut Vec<f64>,
    allele_mapq: &mut Vec<MapqPool>,
    sample_columns: &mut [SampleColumn],
    flank_deletions: &mut Vec<FlankDeletion>,
) {
    // The reference allele is never dropped: it is the bases the record's span spells, and
    // `VcfRecord::new` refuses a record without one.
    let mut kept = vec![false; alleles.len()];
    kept[0] = true;
    let mut any_named = false;
    for column in sample_columns.iter() {
        if let Some(genotype) = column.call.genotype() {
            for allele in genotype.alleles() {
                any_named = true;
                // PANIC-FREE: `VcfRecord::new` asserts the same bound a few lines later, and
                // every genotype here was built over this very allele table.
                kept[usize::from(allele.get())] = true;
            }
        }
    }
    if !any_named || kept.iter().all(|k| *k) {
        return;
    }

    // Old index -> new index, built once and read by all four parallel arrays.
    let mut renumbered = vec![None; alleles.len()];
    let mut next = 0_usize;
    for (old, keep) in kept.iter().enumerate() {
        if *keep {
            renumbered[old] = Some(next);
            next += 1;
        }
    }

    for column in sample_columns.iter_mut() {
        let mut reads = Vec::with_capacity(next);
        let mut dropped_reads = 0_u32;
        for (old, count) in column.read_counts.allele_reads().iter().enumerate() {
            if kept[old] {
                reads.push(*count);
            } else {
                dropped_reads = dropped_reads.saturating_add(*count);
            }
        }
        let unexplained = column
            .read_counts
            .unexplained_reads()
            .saturating_add(dropped_reads);
        column.read_counts = SampleReadCounts::new(reads, unexplained);
        if let SampleCall::Called { genotype, .. } = &mut column.call {
            let renamed: Vec<AlleleId> = genotype
                .alleles()
                .iter()
                .map(|allele| {
                    // PANIC-FREE: every allele a genotype names was marked kept above, so
                    // its slot in `renumbered` is `Some`.
                    let new = renumbered[usize::from(allele.get())]
                        .expect("a genotype's allele is kept by construction");
                    AlleleId(new as u16)
                })
                .collect();
            *genotype = Genotype::new(renamed);
        }
    }

    let mut old = 0_usize;
    alleles.retain(|_| {
        let keep = kept[old];
        old += 1;
        keep
    });
    let mut old = 0_usize;
    expected_copies.retain(|_| {
        let keep = kept[old];
        old += 1;
        keep
    });
    let mut old = 0_usize;
    allele_mapq.retain(|_| {
        let keep = kept[old];
        old += 1;
        keep
    });
    let mut old = 0_usize;
    flank_deletions.retain(|_| {
        let keep = kept[old];
        old += 1;
        keep
    });
}

fn written_call(call: &SampleGenotypeCall) -> SampleCall {
    match call {
        SampleGenotypeCall::Missing => SampleCall::NoCall,
        SampleGenotypeCall::Called {
            genotype,
            genotype_quality,
            reads_were_uninformative,
        } => {
            if *reads_were_uninformative {
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
