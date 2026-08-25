//! The STR row — one sample's genotype likelihoods at one repeat tract.
//!
//! *How well does each candidate genotype explain what this sample's reads showed at this
//! tract?* The answer is one number per genotype, and this module computes the whole row.
//!
//! Everything about **how a single read is scored against a single allele** lives behind the
//! emission seam ([`super::ssr_emission`]); everything about **how those scores are combined
//! into a genotype's likelihood** lives here. That split is what lets the model comparison
//! behind spec §4.1 swap one model for another without touching the arithmetic around it.
//!
//! # The formula, in words
//!
//! A read arrived one of two ways: it was copied from one of this individual's own copies of
//! the tract, or it did not come from this individual's copies at all. Writing that for a
//! genotype `g` whose copy counts are `k_a` over a ploidy `P` (spec §2.1, §4.5):
//!
//! ```text
//! log Lg(g)  =  Σ_o  n_o · log[ (1 − λ) · Σ_a (k_a / P) · Lr(o | a)  +  λ · U ]
//! ```
//!
//! - `o` runs over the sample's observations at this tract and `n_o` is how many reads showed
//!   each one;
//! - `k_a / P` is the chance a read was copied from a copy carrying allele `a`;
//! - `Lr(o | a)` is the emission — [`SsrEmissionModel::emission`] for a read that spanned the
//!   tract, [`SsrEmissionModel::censored_emission`] for one that ran off its own end;
//! - `λ` is the chance the read came from somewhere else entirely, and `U` is what such a read
//!   shows: uniform over the tract lengths the model's support can reach.
//!
//! **The second term is what stops one strange read from ruling out every genotype.** A read
//! from a paralogous tract, a chimera, a somatic length in a long tract — its emission is zero
//! under every candidate, so the bracket collapses to `λ · U` for every genotype alike and the
//! term contributes the same number to every entry of the row.
//!
//! # What a row costs, and why that is the design
//!
//! **`observations × candidates` emission calls, not `observations × genotypes`.** Each
//! emission is computed once into [`SsrRowScratch`]'s cache and read by every genotype that
//! carries the candidate. At six candidates and a diploid that is 6 calls an observation
//! rather than 21 — and spec §8 calls it the design rather than an optimisation, which is why
//! the cache is a field of the caller's scratch rather than a local here.
//!
//! # What this module does not do
//!
//! **Contamination.** Spec §4.5.1's third term — a read from another individual, which is not
//! junk because it shows a length that is a real allele in somebody — is the next step's, and
//! with it the per-locus computation of how many lengths the outlier weight spreads over.
//! Until then both arrive as parameters.

use super::ssr_emission::{SsrCandidate, SsrEmissionModel, SsrScoringContext};
use super::{SsrRowScratch, SsrSampleEvidence};
use crate::ng::calling::GenotypeTableView;
use crate::ng::types::{LogProb, ReadGroupId};

/// How often a read at a repeat tract came from somewhere other than this individual's copies
/// of it — **inherited from production at 0.01 and declared inherited, not fitted.**
///
/// Production sets it here ([`em.rs`](../../../../src/ssr/cohort/em.rs)) and ng keeps the
/// number. It has no source in the parameters fit — nothing measures it — and spec §1.2
/// requires that be said rather than left blank, so this is a named constant awaiting a
/// measurement rather than a finding.
pub const DEFAULT_OUTLIER_WEIGHT: f64 = 0.01;

/// The scoring parameters for every `(read group, candidate)` pair at one locus, as one
/// checked table.
///
/// **The stride has one spelling, which is the whole point of the type.** A context is built
/// per `(read group, candidate)` — a read's chance of slipping is a property of the tract it
/// was copied from, so a 6-repeat candidate and a 12-repeat one at the same locus are drawn
/// from different strata (spec §4.4) — and a bare slice with the indexing written at the call
/// site is two chances to write it differently. Getting it wrong reads a real context that
/// belongs to another candidate: a plausible number, no panic.
#[derive(Debug, Clone, Copy)]
pub struct SsrScoringContextTable<'a> {
    contexts: &'a [SsrScoringContext<'a>],
    candidates: usize,
}

impl<'a> SsrScoringContextTable<'a> {
    /// Wrap the table, checking it is rectangular.
    ///
    /// # Panics
    ///
    /// If the slice is not `read_groups × candidates` entries. A short table would otherwise
    /// surface at whichever locus first reached past its end, or never.
    #[must_use]
    pub fn new(contexts: &'a [SsrScoringContext<'a>], candidates: usize) -> Self {
        assert!(
            candidates > 0,
            "a locus is called over at least its reference allele"
        );
        assert!(
            contexts.len().is_multiple_of(candidates),
            "a context table holds one entry per (read group, candidate): {} entries is not a \
             whole number of rows of {candidates}",
            contexts.len()
        );
        Self {
            contexts,
            candidates,
        }
    }

    /// How many read groups this table covers.
    #[must_use]
    pub fn read_group_count(&self) -> usize {
        self.contexts.len() / self.candidates
    }

    /// How many candidates each read group's row covers.
    #[must_use]
    pub fn candidate_count(&self) -> usize {
        self.candidates
    }

    /// The context for one `(read group, candidate)`.
    ///
    /// # Panics
    ///
    /// **In release as well as debug**, on a read group or candidate past what the table
    /// covers — because the alternative is scoring a read against another candidate's stutter
    /// parameters and getting a number back.
    #[must_use]
    pub fn of(&self, read_group: ReadGroupId, candidate: usize) -> &SsrScoringContext<'a> {
        let group = read_group.get() as usize;
        assert!(
            candidate < self.candidates,
            "candidate {candidate} is past the {} this table covers",
            self.candidates
        );
        assert!(
            group < self.read_group_count(),
            "read group {group} is past the {} this table covers",
            self.read_group_count()
        );
        &self.contexts[group * self.candidates + candidate]
    }
}

/// What the **locus** contributes to a row — the same for every sample called at this tract,
/// and none of it derivable from the genotype table.
///
/// **Grouped rather than passed loose**, and not only to shorten the signature: these four
/// travel together, they are built once per locus by the calling loop, and the next step adds a
/// fifth to them (spec §4.5.1's contamination). A caller that has one of them has all of them.
#[derive(Debug, Clone, Copy)]
pub struct SsrLocusParameters<'a> {
    /// The candidates, in genotype-table allele order.
    ///
    /// **Built [`SsrCandidate`]s rather than bases, because a candidate's repeat count is not
    /// derivable from its bases**: an interrupted tract's byte length divided by the period is
    /// not how many repeats it holds. The locus generator has already measured it, and
    /// re-measuring here would be the duplication spec §7 puts on the alignment module's side
    /// of the boundary.
    pub candidates: &'a [SsrCandidate<'a>],
    /// The stutter and substitution parameters, per `(read group, candidate)` — never hoisted
    /// out of the candidate loop (spec §4.4).
    pub contexts: SsrScoringContextTable<'a>,
    /// How often a read came from somewhere other than this individual's copies of the tract.
    /// [`DEFAULT_OUTLIER_WEIGHT`] is the value to pass.
    pub outlier_weight: f64,
    /// How many tract lengths the outlier weight is spread over — **a property of the candidate
    /// set and the two cutoffs, with no cohort in it** (spec §4.5). A parameter until the next
    /// step computes it here.
    pub reachable_length_count: u32,
}

/// **One sample's log-likelihood for every candidate genotype at one repeat tract**, written
/// into `out` in genotype-table order.
///
/// The module's own documentation carries the formula and what each term is for. This is what
/// a caller has to get right.
///
/// What the locus contributes is [`SsrLocusParameters`], which carries its own documentation for
/// each field and why none of them is derived here.
///
/// # Panics
///
/// On any mismatch between the tables handed in — a row of the wrong width, a context table
/// for a different locus, a ploidy past [`MAX_PLOIDY_COPIES`], or a
/// `reachable_length_count` of zero. **All of them are checked once here rather than per
/// observation**, so that the accessors inside the loops cannot fail in a run; a lazily
/// checked mismatch would surface at whichever locus first reached past an end, or never.
pub fn genotype_log_likelihood_row<Model: SsrEmissionModel>(
    model: &Model,
    evidence: &SsrSampleEvidence<'_>,
    locus: SsrLocusParameters<'_>,
    genotypes: &GenotypeTableView<'_>,
    out: &mut [LogProb],
    scratch: &mut SsrRowScratch<Model::Scratch>,
) {
    let SsrLocusParameters {
        candidates,
        contexts,
        outlier_weight,
        reachable_length_count,
    } = locus;
    let genotype_count = genotypes.genotype_count();
    let allele_count = genotypes.allele_count();
    assert_eq!(
        out.len(),
        genotype_count,
        "a row holds one entry per candidate genotype — {genotype_count}, not {}",
        out.len()
    );
    assert_eq!(
        candidates.len(),
        allele_count,
        "this locus was handed {} candidates and a genotype table over {allele_count} alleles, \
         so one of them belongs to a different locus",
        candidates.len()
    );
    assert_eq!(
        contexts.candidate_count(),
        allele_count,
        "the context table covers {} candidates and this locus is called over {allele_count}",
        contexts.candidate_count()
    );
    assert!(
        reachable_length_count > 0,
        "the outlier term is spread over the lengths the model can reach, and there is always \
         at least one — a candidate's own"
    );
    // **The outlier weight is a share of the reads, so it lives strictly inside 0 and 1**, and
    // both ends are checked because both produce a number rather than a crash. Above one,
    // `1 − λ` is negative and the bracket goes negative wherever a genotype explains the read,
    // so `ln` returns `NaN` — measured at λ = 1.5, two of three genotypes came back `NaN` and
    // the third a plausible −12.74. At exactly zero the row loses its only floor and a read
    // nothing explains takes every genotype to `−∞`, whose differences are `NaN` — which is
    // precisely what spec §4.5's junk term exists to prevent. A measurement that genuinely
    // wants no outlier term should say so rather than reach it through this argument.
    assert!(
        outlier_weight > 0.0 && outlier_weight < 1.0,
        "the outlier weight is the share of reads that came from somewhere else, so it lies \
         strictly inside 0 and 1 — not {outlier_weight}"
    );

    // **The cache is filled with `NaN` and not with zero.** An unwritten slot has to hold
    // something the row cannot mistake for a real score, and zero is exactly that mistake: a
    // slip a candidate cannot reach legitimately scores zero (spec §4.2), so zeros would make
    // *never computed* and *computed as impossible* the same value.
    scratch.prepare_emissions(evidence, allele_count, f64::NAN);
    fill_emissions(model, evidence, candidates, contexts, scratch);

    // `k / P` for every copy count a genotype can carry, shared with the SNP/indel row so the
    // two cannot disagree about what that is. It also carries the ploidy check.
    let copy_share = super::copy_shares(genotypes.ploidy());

    // The outlier half of the mixture: a read from somewhere else shows any reachable length,
    // and one length is as likely as another (spec §4.5). It is the same number for every
    // observation and every genotype, so it is computed once.
    let from_the_junk_distribution = outlier_weight / f64::from(reachable_length_count);
    let from_this_individual = 1.0 - outlier_weight;

    for slot in out.iter_mut() {
        *slot = LogProb(0.0);
    }

    // **The observations are walked in the order the caller handed them**, and this row
    // imposes none of its own — no sorting, no bucketing, no re-grouping. The caller always
    // hands it the merge's order, and that is what makes a run reproducible at any worker
    // count (spec §12 test 8).
    let counts = genotypes.genotype_allele_counts();
    for (position, observation) in evidence.observations.iter().enumerate() {
        let reads = f64::from(observation.num_obs);
        for (genotype, slot) in out.iter_mut().enumerate() {
            let carried_copies = &counts[genotype * allele_count..][..allele_count];
            // **The copy-weighted mixture over the genotype's own alleles.** A candidate no
            // copy carries is skipped rather than multiplied by zero, which is what keeps the
            // cost proportional to the ploidy rather than to the candidate count.
            let mut explained_by_this_genotype = 0.0;
            for (candidate, &copies) in carried_copies.iter().enumerate() {
                if copies == 0 {
                    continue;
                }
                explained_by_this_genotype +=
                    copy_share[copies as usize] * scratch.emission_at(position, candidate);
            }
            slot.0 += reads
                * (from_this_individual * explained_by_this_genotype + from_the_junk_distribution)
                    .ln();
        }
    }
}

/// Score every `(observation, candidate)` pair into the cache, routing each observation by
/// what its reads actually witnessed.
///
/// **The witness decides the method, and nothing else does.** A read that spanned the tract
/// pins a length and goes to [`SsrEmissionModel::emission`]; a read that ran off its own end
/// proves only a lower bound and goes to [`SsrEmissionModel::censored_emission`]. Scoring a
/// partial as though it were complete mis-scores it as a *short* allele, because its bases are
/// a prefix of the truth (spec §5.1) — so the split is taken from
/// [`SsrSampleEvidence`]'s two filters, which are the single place that decides what reaches
/// the censored term, rather than re-derived from the bases here.
///
/// **The two filters together cover every observation exactly once**, which is what lets the
/// cache be filled with `NaN` and still hold no `NaN` when this returns: they enumerate the
/// same slice and split it on an exhaustive match.
fn fill_emissions<Model: SsrEmissionModel>(
    model: &Model,
    evidence: &SsrSampleEvidence<'_>,
    candidates: &[SsrCandidate<'_>],
    contexts: SsrScoringContextTable<'_>,
    scratch: &mut SsrRowScratch<Model::Scratch>,
) {
    for (position, observation) in evidence.complete_observations() {
        for (candidate, allele) in candidates.iter().enumerate() {
            let context = contexts.of(observation.read_group, candidate);
            let scored = model.emission(
                &observation.bases,
                allele,
                context,
                scratch.model_scratch_mut(),
            );
            scratch.set_emission(position, candidate, scored);
        }
    }

    for (position, observation) in evidence.partial_observations() {
        for (candidate, allele) in candidates.iter().enumerate() {
            let context = contexts.of(observation.read_group, candidate);
            let scored = model.censored_emission(
                &observation.bases,
                allele,
                context,
                scratch.model_scratch_mut(),
            );
            scratch.set_emission(position, candidate, scored);
        }
    }
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;
    use std::num::NonZeroU32;

    use super::*;
    use crate::ng::alignment::StutterModel;
    use crate::ng::calling::genotype_table::GenotypeTable;
    use crate::ng::calling::likelihood::ssr_emission::{
        StutterSubstitutionEmission, StutterSubstitutionScratch,
    };
    use crate::ng::calling::likelihood::stutter_rates::stutter_model_for;
    use crate::ng::locus_generation::{LocusLen, ReadWitness, SequenceObservation, SsrDetail};
    use crate::ng::parameter_estimation::Provenance;
    use crate::ng::parameter_estimation::joint::ssr_fit::Slippage;
    use crate::ng::types::{ErrorRate, Motif, Ploidy};

    /// The per-base substitution rate the fixtures score at — **not** a floating-point
    /// tolerance, which is what `EPSILON` would read as a few lines from `f64::EPSILON`.
    const SUBSTITUTION_RATE: f64 = 1e-3;

    fn a_motif(bases: &[u8]) -> Motif {
        Motif::new(bases).expect("a valid test motif")
    }

    fn repeats(count: u32) -> NonZeroU32 {
        NonZeroU32::new(count).expect("a test candidate always holds a repeat")
    }

    /// A slippage row at a stated level — what makes two candidates' contexts different.
    fn a_model_slipping(level: f64) -> StutterModel {
        stutter_model_for(&Slippage {
            level,
            shorter_share: 0.83,
            fall_off: 0.35,
        })
    }

    fn a_tract(motif: &[u8], repeat_count: usize) -> Vec<u8> {
        motif.repeat(repeat_count)
    }

    /// One observation: some bases, seen by `reads` reads that spanned the whole tract.
    fn spanning(bases: &[u8], reads: u32) -> SequenceObservation {
        SequenceObservation {
            bases: bases.to_vec().into_boxed_slice(),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(0),
            num_obs: reads,
            num_fwd: reads,
            q_sum: -10.0 * f64::from(reads),
            mapq_sum: 60 * reads,
            mapq_sum_sq: u64::from(reads) * 3_600,
            placed_left: 0,
            chain_ids: Vec::new(),
        }
    }

    /// The detail an STR locus carries — the motif, and flanks nothing here reads.
    fn a_detail(motif: &[u8]) -> SsrDetail {
        SsrDetail {
            motif: a_motif(motif),
            left_flank: Box::from(&b"GGGGG"[..]),
            right_flank: Box::from(&b"TTTTT"[..]),
        }
    }

    /// Everything a row needs at one locus, owned so the borrows below stay simple.
    ///
    /// **Every candidate gets its own [`StutterModel`] and its own substitution rate**, and
    /// that is not decoration. A read's chance of slipping is a property of the tract it was
    /// copied from, so candidates of different repeat counts are drawn from different strata
    /// and the row must look each one up separately (spec §4.4). A fixture that gave every
    /// candidate the same parameters could not tell that apart from a row that read candidate
    /// zero's context for all of them — and the first version of this file could not: hoisting
    /// the lookup out of the candidate loop left every row bit-identical.
    ///
    /// **The read groups differ too**, for the same reason on the other axis of the table.
    struct Fixture {
        motif: Motif,
        /// One per `(read group, candidate)`, in the table's own order.
        models: Vec<StutterModel>,
        /// One per `(read group, candidate)`, matching `models`.
        substitution_rates: Vec<f64>,
        read_groups: usize,
        candidate_bases: Vec<Vec<u8>>,
        repeat_counts: Vec<u32>,
    }

    impl Fixture {
        /// Candidates of `repeat_counts` whole copies of `motif`, in that order, scored by one
        /// read group.
        fn of(motif: &[u8], repeat_counts: &[u32]) -> Self {
            Self::of_groups(motif, repeat_counts, 1)
        }

        /// The same, across `read_groups` read groups — each with its own slippage row, as a
        /// per-chemistry fit gives them.
        fn of_groups(motif: &[u8], repeat_counts: &[u32], read_groups: usize) -> Self {
            let mut models = Vec::new();
            let mut substitution_rates = Vec::new();
            for group in 0..read_groups {
                for count in repeat_counts {
                    // Longer tracts slip more, and one lane slips more than another: two
                    // separate reasons for two contexts to differ, so a row that dropped
                    // either axis of the lookup has something to fail against.
                    //
                    // **Both are keyed by the candidate's repeat count, never by its position
                    // in the table** — that is what a stratum lookup does, and a fixture keyed
                    // by position would make permuting the candidates change the answer for a
                    // reason that has nothing to do with the row.
                    models.push(a_model_slipping(
                        0.01 + 0.004 * f64::from(*count) + 0.002 * group as f64,
                    ));
                    substitution_rates
                        .push(SUBSTITUTION_RATE * (1.0 + 0.1 * f64::from(*count) + group as f64));
                }
            }
            Self {
                motif: a_motif(motif),
                models,
                substitution_rates,
                read_groups,
                candidate_bases: repeat_counts
                    .iter()
                    .map(|count| a_tract(motif, *count as usize))
                    .collect(),
                repeat_counts: repeat_counts.to_vec(),
            }
        }

        fn candidates(&self) -> Vec<SsrCandidate<'_>> {
            self.candidate_bases
                .iter()
                .zip(&self.repeat_counts)
                .map(|(bases, count)| SsrCandidate {
                    bases,
                    repeat_count: repeats(*count),
                })
                .collect()
        }

        /// The whole `(read group, candidate)` table, in the order
        /// [`SsrScoringContextTable`] indexes it.
        fn contexts<'a>(&'a self, candidates: &[SsrCandidate<'_>]) -> Vec<SsrScoringContext<'a>> {
            let mut contexts = Vec::new();
            for group in 0..self.read_groups {
                for (candidate, allele) in candidates.iter().enumerate() {
                    let at = group * candidates.len() + candidate;
                    contexts.push(SsrScoringContext::new(
                        &self.motif,
                        &self.models[at],
                        allele,
                        ErrorRate::try_new(self.substitution_rates[at]).expect("a valid rate"),
                        [Provenance::FittedHere],
                    ));
                }
            }
            contexts
        }
    }

    /// Build the row for one set of observations, at a stated ploidy and outlier weight.
    fn score_row_at(
        fixture: &Fixture,
        observations: &[SequenceObservation],
        ploidy: u8,
        reachable_lengths: u32,
        outlier_weight: f64,
    ) -> Vec<LogProb> {
        let detail = a_detail(fixture.motif.as_bytes());
        let evidence = SsrSampleEvidence::new(observations, &detail);
        let candidates = fixture.candidates();
        let contexts = fixture.contexts(&candidates);
        let table = GenotypeTable::build(
            Ploidy::try_new(ploidy).expect("a valid ploidy"),
            candidates.len(),
        );
        let view = table.view();
        let mut out = vec![LogProb(f64::NAN); view.genotype_count()];
        let mut scratch = SsrRowScratch::<StutterSubstitutionScratch>::default();
        genotype_log_likelihood_row(
            &StutterSubstitutionEmission,
            &evidence,
            SsrLocusParameters {
                candidates: &candidates,
                contexts: SsrScoringContextTable::new(&contexts, candidates.len()),
                outlier_weight,
                reachable_length_count: reachable_lengths,
            },
            &view,
            &mut out,
            &mut scratch,
        );
        out
    }

    /// The same at the inherited outlier weight, which is what every test but two wants.
    fn score_row(
        fixture: &Fixture,
        observations: &[SequenceObservation],
        ploidy: u8,
        reachable_lengths: u32,
    ) -> Vec<LogProb> {
        score_row_at(
            fixture,
            observations,
            ploidy,
            reachable_lengths,
            DEFAULT_OUTLIER_WEIGHT,
        )
    }

    /// Whether no candidate at all can produce this observation — what makes a read *junk*
    /// rather than weak evidence, and the thing the cancellation test is about.
    fn every_candidate_scores_it_zero(
        fixture: &Fixture,
        observation: &SequenceObservation,
    ) -> bool {
        let candidates = fixture.candidates();
        let contexts = fixture.contexts(&candidates);
        let mut scratch = StutterSubstitutionScratch::default();
        candidates
            .iter()
            .zip(&contexts)
            .all(|(candidate, context)| {
                StutterSubstitutionEmission.emission(
                    &observation.bases,
                    candidate,
                    context,
                    &mut scratch,
                ) == 0.0
            })
    }

    /// **Spec §12's sixth test: the junk term cancels for a read nothing explains.**
    ///
    /// A read from a paralogous tract, a chimera, a somatic length — its emission is zero under
    /// every candidate, so the bracket collapses to `λ · U`, the same number under every
    /// genotype. **What that must not do is change how far apart two genotypes are**, because
    /// that is what the caller normalises and calls with.
    ///
    /// # It cancels to one unit in the last place of the entries, not of their difference
    ///
    /// The specification asked for bitwise, and no implementation of this formula can give it:
    /// the junk term is added to each genotype's running total, and `(a + k) − (b + k)` is not
    /// `a − b` in floating point however carefully `k` is computed.
    ///
    /// **The unit matters more than the number here, and getting it wrong is how this test was
    /// first written.** Counting units in the last place *of the difference* measures the true
    /// rounding error scaled by `|entry| / |separation|`, and that ratio is set by how many junk
    /// reads there are — so the same fixture reports 16 units at 3 junk reads and 3,072 at 300,
    /// with nothing about the row having changed. Measured relative to **the entries' own
    /// magnitude**, which is where the rounding actually happens, the worst disagreement over
    /// every cell below is one `f64::EPSILON` and stays there. That is the same shape of bound
    /// `permuting_the_observations_and_the_candidates_moves_no_genotype_meaningfully` uses, and
    /// what spec §12's eighth test means by "the same relative bound as test 9".
    ///
    /// The sweep varies the junk read count for exactly that reason: a fixture with three junk
    /// reads cannot tell a robust bound from a lucky one.
    #[test]
    fn a_read_nothing_explains_moves_no_genotype_against_another() {
        let mut worst_relative = 0.0f64;
        let mut cells = 0usize;
        let mut cells_that_moved = 0usize;

        for candidate_counts in [vec![4u32, 5], vec![4, 5, 6, 7], vec![3, 6, 9]] {
            let fixture = Fixture::of(b"CA", &candidate_counts);
            let longest = *candidate_counts.iter().max().expect("a candidate");
            let explained = spanning(&a_tract(b"CA", candidate_counts[0] as usize), 7);
            let other = spanning(&a_tract(b"CA", longest as usize), 4);

            for junk_reads in [3u32, 30, 100, 300] {
                // **Eleven repeats past the *longest* candidate, not the shortest.** One past
                // the whole-repeat cutoff of a short candidate is still within the cutoff of a
                // longer one, which the distribution scores rather than refuses — so a junk
                // read built from the shortest candidate is not junk at all. The first draft of
                // this fixture made that mistake and failed for a reason that had nothing to do
                // with the property under test.
                let junk = spanning(&a_tract(b"CA", (longest + 11) as usize), junk_reads);

                // **The fixture has to be junk, or this test is about nothing.** Asserted
                // rather than reasoned about, because the reasoning is what went wrong.
                assert!(
                    every_candidate_scores_it_zero(&fixture, &junk),
                    "some candidate explains the junk read, so the term under test does not \
                     cancel"
                );

                let without = score_row(&fixture, &[explained.clone(), other.clone()], 2, 12);
                let with = score_row(&fixture, &[explained.clone(), other.clone(), junk], 2, 12);

                assert!(
                    with.iter().zip(&without).all(|(a, b)| a.0 < b.0),
                    "the junk read should cost every genotype something"
                );
                for first in 0..without.len() {
                    for second in 0..without.len() {
                        let moved = (with[first].0 - with[second].0)
                            - (without[first].0 - without[second].0);
                        // The entries are where the rounding happens, so they are what the
                        // error is measured against.
                        let magnitude = with[first].0.abs().max(with[second].0.abs());
                        let relative = moved.abs() / magnitude;
                        assert!(
                            relative <= f64::EPSILON,
                            "genotypes {first} and {second} moved {moved} at {junk_reads} junk \
                             reads — {relative} of an entry of magnitude {magnitude}"
                        );
                        worst_relative = worst_relative.max(relative);
                        if moved != 0.0 {
                            cells_that_moved += 1;
                        }
                        cells += 1;
                    }
                }
            }
        }

        // **The sweep has to reach cells where the two genuinely differ**, or the bound is
        // being read off arithmetic that happened to cancel and this test would pass just as
        // well against the bitwise claim it replaces.
        assert_eq!(cells, 4 * (9 + 100 + 36), "the sweep changed size");
        assert!(
            cells_that_moved > 0,
            "no cell moved at all, so this sweep cannot tell a real bound from a bitwise one"
        );
        assert!(
            worst_relative > 0.0 && worst_relative <= f64::EPSILON,
            "the worst relative disagreement was {worst_relative}"
        );
    }

    /// **Spec §12's seventh test: ploidy generality, and both of its cases.**
    ///
    /// At ploidy 4 with every read matching one allele, a genotype carrying two copies each of
    /// two alleles scores **between** the two homozygous quadruples. Where the reads are split
    /// between those alleles it scores **above both** — and that second case is the whole reason
    /// a mixed genotype is callable, so pinning only the first would leave a wrong copy
    /// weighting undetected.
    #[test]
    fn a_mixed_genotype_sits_between_the_homozygotes_and_above_them_when_the_reads_split() {
        let fixture = Fixture::of(b"CA", &[4, 5]);
        let four = a_tract(b"CA", 4);
        let five = a_tract(b"CA", 5);

        let table = GenotypeTable::build(Ploidy::try_new(4).expect("a valid ploidy"), 2);
        let view = table.view();
        let counts = view.genotype_allele_counts();
        let find = |first: u32, second: u32| {
            (0..view.genotype_count())
                .find(|genotype| {
                    counts[genotype * 2] == first && counts[genotype * 2 + 1] == second
                })
                .expect("the genotype table holds every copy split")
        };
        let all_four = find(4, 0);
        let all_five = find(0, 4);
        let mixed = find(2, 2);

        // Every read matches the four-repeat allele.
        let matching = score_row(&fixture, &[spanning(&four, 10)], 4, 12);
        assert!(
            matching[mixed].0 < matching[all_four].0 && matching[mixed].0 > matching[all_five].0,
            "with every read on one allele the mixed genotype must sit between the two \
             homozygous quadruples: {:?}",
            (
                matching[all_four].0,
                matching[mixed].0,
                matching[all_five].0
            )
        );

        // The reads split evenly between the two alleles.
        let split = score_row(&fixture, &[spanning(&four, 5), spanning(&five, 5)], 4, 12);
        assert!(
            split[mixed].0 > split[all_four].0 && split[mixed].0 > split[all_five].0,
            "with the reads split the mixed genotype must beat both homozygous quadruples: {:?}",
            (split[all_four].0, split[mixed].0, split[all_five].0)
        );
    }

    /// **Spec §12's eighth test: the row imposes no order of its own.**
    ///
    /// It must not sort, bucket or re-group behind the caller: the caller always hands it the
    /// merge's order, and that is what makes a run reproducible at any worker count.
    ///
    /// **Not asserted bitwise, and the specification says why**: permuting the observations
    /// *is* changing the summation order, so the two rows may differ in the last bits. What is
    /// asserted is that they differ by no more than that — measured here at zero units in the
    /// last place, which is stronger than the bound and is recorded rather than relied on.
    #[test]
    fn permuting_the_observations_and_the_candidates_moves_no_genotype_meaningfully() {
        let fixture = Fixture::of(b"CA", &[4, 5, 6]);
        let observations = [
            spanning(&a_tract(b"CA", 4), 7),
            spanning(&a_tract(b"CA", 5), 3),
            spanning(&a_tract(b"CA", 6), 2),
        ];
        let forwards = score_row(&fixture, &observations, 2, 12);

        let mut backwards_observations = observations.to_vec();
        backwards_observations.reverse();
        let backwards = score_row(&fixture, &backwards_observations, 2, 12);

        for (genotype, (one, other)) in forwards.iter().zip(&backwards).enumerate() {
            let apart = (one.0 - other.0).abs();
            assert!(
                apart <= 2.0 * f64::EPSILON * one.0.abs().max(1.0),
                "genotype {genotype} moved {apart} nats when the observations were permuted"
            );
        }

        // **And the candidates, which the first version of this test did not do despite its
        // name.** Reversing the candidate order relabels the genotypes, so the match is made on
        // copy counts rather than by mirroring the index — at three candidates the reversal
        // maps genotype indices 0→5, 1→4, 2→2, 3→3, 4→1, 5→0, and a naive mirror compares two
        // genotypes that are twenty nats apart.
        let reversed_counts: Vec<u32> = vec![6, 5, 4];
        let reversed_fixture = Fixture::of(b"CA", &reversed_counts);
        let reversed = score_row(&reversed_fixture, &observations, 2, 12);

        let table = GenotypeTable::build(Ploidy::try_new(2).expect("a valid ploidy"), 3);
        let view = table.view();
        let counts = view.genotype_allele_counts();
        let mut matched = 0usize;
        for forward_genotype in 0..view.genotype_count() {
            let carried = &counts[forward_genotype * 3..][..3];
            // The same genotype under the reversed candidate order carries the mirrored counts.
            let mirrored: Vec<u32> = carried.iter().rev().copied().collect();
            let reversed_genotype = (0..view.genotype_count())
                .find(|genotype| counts[genotype * 3..][..3] == mirrored[..])
                .expect("every copy split appears under either ordering");
            assert_eq!(
                forwards[forward_genotype].0.to_bits(),
                reversed[reversed_genotype].0.to_bits(),
                "genotype {carried:?} moved when the candidates were permuted"
            );
            matched += 1;
        }
        assert_eq!(matched, 6, "the genotype table changed shape");
    }

    /// **A row costs `observations × candidates` emission calls — not `× genotypes`.**
    ///
    /// Spec §8 calls that the design rather than an optimisation, and this is the only test
    /// that can tell the difference: every other one would pass just as well if the row
    /// recomputed each emission for every genotype, because the numbers would be identical.
    ///
    /// The count is instrumented rather than argued, and the fixture is chosen so the three
    /// plausible costs are three different numbers. At three observations, three candidates and
    /// a diploid — six genotypes, nine carried-allele slots — **the design costs 9 calls**,
    /// recomputing per genotype would cost 18, and recomputing per carried allele would cost
    /// 27. An earlier version of this comment gave 18 the "per carried allele" rule, which is
    /// the arithmetic for a different shape.
    #[test]
    fn one_row_scores_each_observation_against_each_candidate_exactly_once() {
        /// Counts what the row asks of it and forwards to the real model.
        struct Counting {
            inner: StutterSubstitutionEmission,
            complete: Cell<usize>,
            censored: Cell<usize>,
        }

        impl SsrEmissionModel for Counting {
            type Scratch = StutterSubstitutionScratch;

            fn emission(
                &self,
                observation: &[u8],
                candidate: &SsrCandidate<'_>,
                context: &SsrScoringContext<'_>,
                scratch: &mut Self::Scratch,
            ) -> f64 {
                self.complete.set(self.complete.get() + 1);
                self.inner
                    .emission(observation, candidate, context, scratch)
            }

            fn censored_emission(
                &self,
                witnessed_prefix: &[u8],
                candidate: &SsrCandidate<'_>,
                context: &SsrScoringContext<'_>,
                scratch: &mut Self::Scratch,
            ) -> f64 {
                self.censored.set(self.censored.get() + 1);
                self.inner
                    .censored_emission(witnessed_prefix, candidate, context, scratch)
            }
        }

        let fixture = Fixture::of(b"CA", &[4, 5, 6]);
        let observations = [
            spanning(&a_tract(b"CA", 4), 7),
            spanning(&a_tract(b"CA", 5), 3),
            spanning(&a_tract(b"CA", 6), 2),
        ];
        let detail = a_detail(b"CA");
        let evidence = SsrSampleEvidence::new(&observations, &detail);
        let candidates = fixture.candidates();
        let contexts = fixture.contexts(&candidates);
        let table = GenotypeTable::build(Ploidy::try_new(2).expect("a valid ploidy"), 3);
        let view = table.view();
        assert_eq!(
            view.genotype_count(),
            6,
            "the fixture must hold more genotypes than candidates, or it cannot tell the two \
             costs apart"
        );

        let model = Counting {
            inner: StutterSubstitutionEmission,
            complete: Cell::new(0),
            censored: Cell::new(0),
        };
        let mut out = vec![LogProb(f64::NAN); view.genotype_count()];
        let mut scratch = SsrRowScratch::<StutterSubstitutionScratch>::default();
        genotype_log_likelihood_row(
            &model,
            &evidence,
            SsrLocusParameters {
                candidates: &candidates,
                contexts: SsrScoringContextTable::new(&contexts, candidates.len()),
                outlier_weight: DEFAULT_OUTLIER_WEIGHT,
                reachable_length_count: 12,
            },
            &view,
            &mut out,
            &mut scratch,
        );

        assert_eq!(
            model.complete.get(),
            observations.len() * candidates.len(),
            "the row scored each (observation, candidate) pair {} times",
            model.complete.get()
        );
        assert_eq!(model.censored.get(), 0, "no observation here is censored");
    }

    /// **Spec §12's seventh test, first half: a biallelic diploid reproduced by hand.**
    ///
    /// The other ploidy test pins *orderings* — which genotype beats which — and every one of
    /// those survives a copy weighting that is wrong by a constant factor. This one recomputes
    /// the whole formula outside the row, from the two emissions and the seven numbers spec
    /// §2.1 names, and requires the answer to the bit.
    ///
    /// **It is the test that pins `k / P` itself.** Without it, deleting the division by the
    /// ploidy — so the weights sum to the ploidy instead of to one — passes every other test in
    /// this file while moving every entry of the row by about 7 nats, a factor of a thousand in
    /// likelihood.
    #[test]
    fn a_biallelic_diploid_row_matches_the_formula_computed_by_hand() {
        let fixture = Fixture::of(b"CA", &[4, 5]);
        let candidates = fixture.candidates();
        let contexts = fixture.contexts(&candidates);

        let observations = [
            spanning(&a_tract(b"CA", 4), 6),
            spanning(&a_tract(b"CA", 5), 2),
        ];
        // The emissions, taken straight from the model — the row's only input this test does
        // not recompute, because it is the seam's business and not the row's.
        let mut model_scratch = StutterSubstitutionScratch::default();
        let emission = |observation: usize, candidate: usize, scratch: &mut _| -> f64 {
            StutterSubstitutionEmission.emission(
                &observations[observation].bases,
                &candidates[candidate],
                &contexts[candidate],
                scratch,
            )
        };
        let scored: Vec<Vec<f64>> = (0..observations.len())
            .map(|observation| {
                (0..candidates.len())
                    .map(|candidate| emission(observation, candidate, &mut model_scratch))
                    .collect()
            })
            .collect();

        let reachable_lengths = 12u32;
        let junk = DEFAULT_OUTLIER_WEIGHT / f64::from(reachable_lengths);
        let own = 1.0 - DEFAULT_OUTLIER_WEIGHT;
        // The three diploid genotypes over two alleles, in genotype-table order: (2,0), (1,1),
        // (0,2). `k / P` is 1.0, 0.5 and 0.0 for the first allele, and the complement for the
        // second.
        let by_hand: Vec<f64> = [(1.0, 0.0), (0.5, 0.5), (0.0, 1.0)]
            .iter()
            .map(|(first, second)| {
                observations
                    .iter()
                    .enumerate()
                    .map(|(observation, entry)| {
                        let explained =
                            first * scored[observation][0] + second * scored[observation][1];
                        f64::from(entry.num_obs) * (own * explained + junk).ln()
                    })
                    .sum()
            })
            .collect();

        let row = score_row(&fixture, &observations, 2, reachable_lengths);
        assert_eq!(row.len(), by_hand.len(), "the genotype table changed shape");
        for (genotype, (slot, expected)) in row.iter().zip(&by_hand).enumerate() {
            assert_eq!(
                slot.0.to_bits(),
                expected.to_bits(),
                "genotype {genotype}: the row gave {} and the formula {expected}",
                slot.0
            );
        }

        // And the copy weights really are what the fixture claims: a row whose weights summed
        // to the ploidy rather than to one would be about seven nats away from this.
        assert!(
            (row[0].0 - (-11.521_002)).abs() < 1e-5,
            "the fixture moved: {}",
            row[0].0
        );
    }

    /// **A candidate is scored against its own stratum's parameters, not candidate zero's.**
    ///
    /// Spec §4.4 is explicit that the context is built per `(read group, candidate)` and that
    /// the lookup may not be hoisted out of the candidate loop: candidates of different repeat
    /// counts slip at measurably different rates. **Hoisting it is a one-character edit that
    /// nothing else in this file can see** — every other fixture would give the same answer,
    /// because they would all be reading identical parameters.
    ///
    /// The same on the other axis: two read groups with different slippage rows must score the
    /// same bases differently.
    #[test]
    fn each_candidate_and_read_group_is_scored_against_its_own_parameters() {
        let fixture = Fixture::of_groups(b"CA", &[4, 6], 2);
        let candidates = fixture.candidates();
        let contexts = fixture.contexts(&candidates);
        let table = SsrScoringContextTable::new(&contexts, candidates.len());
        let detail = a_detail(b"CA");

        // The table's four entries are four different parameter sets, which is what makes the
        // rows below able to disagree at all.
        assert_eq!(table.read_group_count(), 2);
        for group in 0..2u32 {
            let first = table.of(ReadGroupId(group), 0);
            let second = table.of(ReadGroupId(group), 1);
            assert!(
                first.stutter.same_length_share() != second.stutter.same_length_share(),
                "the two candidates of read group {group} share a slippage row, so this test \
                 cannot see a hoisted lookup"
            );
        }
        assert!(
            table.of(ReadGroupId(0), 0).stutter.same_length_share()
                != table.of(ReadGroupId(1), 0).stutter.same_length_share(),
            "the two read groups share a slippage row"
        );

        // A read of the same bases from the two read groups must be scored differently.
        let bases = a_tract(b"CA", 5);
        let from_first = {
            let mut observation = spanning(&bases, 5);
            observation.read_group = ReadGroupId(0);
            observation
        };
        let from_second = {
            let mut observation = spanning(&bases, 5);
            observation.read_group = ReadGroupId(1);
            observation
        };

        let score_with = |observation: &SequenceObservation| {
            let held = [observation.clone()];
            let evidence = SsrSampleEvidence::new(&held, &detail);
            let genotypes = GenotypeTable::build(
                Ploidy::try_new(2).expect("a valid ploidy"),
                candidates.len(),
            );
            let view = genotypes.view();
            let mut out = vec![LogProb(f64::NAN); view.genotype_count()];
            let mut scratch = SsrRowScratch::<StutterSubstitutionScratch>::default();
            genotype_log_likelihood_row(
                &StutterSubstitutionEmission,
                &evidence,
                SsrLocusParameters {
                    candidates: &candidates,
                    contexts: table,
                    outlier_weight: DEFAULT_OUTLIER_WEIGHT,
                    reachable_length_count: 12,
                },
                &view,
                &mut out,
                &mut scratch,
            );
            out
        };

        let first = score_with(&from_first);
        let second = score_with(&from_second);
        assert!(
            first.iter().zip(&second).any(|(a, b)| a.0 != b.0),
            "the two read groups scored the same read identically, so the row is not reading \
             the read group's own parameters"
        );
    }

    /// **A read that ran out inside the tract reaches the censored term, and is not scored as a
    /// short allele.**
    ///
    /// Without a partial observation in a fixture, three separate defects pass every other test
    /// in this file: routing the partial loop to `emission`, indexing the emission cache by a
    /// dense counter inside the filtered loop, and any change to how a witness is read.
    ///
    /// **The bases are deliberately a prefix of the longer candidate**, which is the case the
    /// whole censored term exists for: scored as a complete read they say *this sample carries
    /// the short allele*, and scored as a lower bound they say *the tract is at least this
    /// long*, which the longer candidate satisfies outright.
    #[test]
    fn a_read_that_ran_out_is_scored_as_a_lower_bound_and_not_as_a_short_allele() {
        let fixture = Fixture::of(b"CA", &[4, 8]);
        let detail = a_detail(b"CA");
        let witnessed = a_tract(b"CA", 4);

        let mut partial = spanning(&witnessed, 9);
        partial.read_witness = ReadWitness::from_left(8, LocusLen::from_positions(16))
            .expect("a partial witness of eight of sixteen positions");
        assert!(
            matches!(partial.read_witness, ReadWitness::Partial { .. }),
            "the fixture must actually be partial"
        );
        let complete = spanning(&witnessed, 9);

        // The cache is indexed by position in the *whole* observation slice, so the partial is
        // put first and a complete observation above it: a dense counter inside the filtered
        // loop addresses the wrong row for that one, and nothing else here would notice.
        let observations = [partial, spanning(&a_tract(b"CA", 8), 4)];
        let evidence = SsrSampleEvidence::new(&observations, &detail);
        assert_eq!(evidence.partial_observations().count(), 1);
        assert_eq!(evidence.complete_observations().count(), 1);

        let as_lower_bound = score_row(&fixture, &observations, 2, 12);
        let as_short_allele = score_row(
            &fixture,
            &[complete, spanning(&a_tract(b"CA", 8), 4)],
            2,
            12,
        );

        // Genotype order over two alleles at ploidy 2: (2,0), (1,1), (0,2).
        let eight_eight = 2;
        let four_four = 0;
        assert!(
            as_lower_bound[eight_eight].0 > as_short_allele[eight_eight].0,
            "read as a lower bound, eight bases of tract should cost the eight-repeat \
             homozygote less than reading it as a whole four-repeat allele: {} against {}",
            as_lower_bound[eight_eight].0,
            as_short_allele[eight_eight].0
        );
        // **What separates the two readings is which allele the read is evidence *for*, not
        // the absolute score.** A lower bound is never less likely than the complete read of the
        // same bases under *any* candidate — that is the censored term's own invariant — so both
        // entries rise. The thing that must move is how far apart the two homozygotes are.
        let separated_as_a_short_allele =
            as_short_allele[four_four].0 - as_short_allele[eight_eight].0;
        let separated_as_a_lower_bound =
            as_lower_bound[four_four].0 - as_lower_bound[eight_eight].0;
        assert!(
            separated_as_a_short_allele > separated_as_a_lower_bound,
            "read as a whole allele the bases should favour the four-repeat homozygote more \
             than the same bases read as a lower bound do: {separated_as_a_short_allele} \
             against {separated_as_a_lower_bound}"
        );
    }

    /// **The row refuses an outlier weight that is not a share of the reads.**
    ///
    /// Both ends produce a number rather than a crash, which is why both are checked: above one
    /// the explained half of the mixture goes negative and `ln` returns `NaN`, and at exactly
    /// zero a read nothing explains takes every genotype to `−∞`, whose differences are `NaN` —
    /// the collapse spec §4.5's junk term exists to prevent.
    #[test]
    #[should_panic(expected = "the outlier weight is the share of reads")]
    fn an_outlier_weight_above_one_is_refused() {
        let fixture = Fixture::of(b"CA", &[4, 5]);
        let observations = [spanning(&a_tract(b"CA", 4), 5)];
        score_row_at(&fixture, &observations, 2, 12, 1.5);
    }

    #[test]
    #[should_panic(expected = "the outlier weight is the share of reads")]
    fn an_outlier_weight_of_zero_is_refused() {
        let fixture = Fixture::of(b"CA", &[4, 5]);
        let observations = [spanning(&a_tract(b"CA", 4), 5)];
        score_row_at(&fixture, &observations, 2, 12, 0.0);
    }

    /// **A locus reaching no lengths at all is refused**, because the outlier weight is spread
    /// over them and there is always at least one — a candidate's own.
    #[test]
    #[should_panic(expected = "the outlier term is spread over the lengths")]
    fn a_locus_reaching_no_lengths_is_refused() {
        let fixture = Fixture::of(b"CA", &[4, 5]);
        let observations = [spanning(&a_tract(b"CA", 4), 5)];
        score_row(&fixture, &observations, 2, 0);
    }

    /// **A candidate set that does not match the genotype table is refused**, rather than
    /// scoring a row over one locus's alleles and another's genotypes.
    #[test]
    #[should_panic(expected = "belongs to a different locus")]
    fn a_candidate_set_the_genotype_table_does_not_cover_is_refused() {
        let fixture = Fixture::of(b"CA", &[4, 5, 6]);
        let detail = a_detail(b"CA");
        let observations = [spanning(&a_tract(b"CA", 4), 5)];
        let evidence = SsrSampleEvidence::new(&observations, &detail);
        let candidates = fixture.candidates();
        let contexts = fixture.contexts(&candidates);
        // A genotype table over two alleles, against three candidates.
        let table = GenotypeTable::build(Ploidy::try_new(2).expect("a valid ploidy"), 2);
        let view = table.view();
        let mut out = vec![LogProb(f64::NAN); view.genotype_count()];
        let mut scratch = SsrRowScratch::<StutterSubstitutionScratch>::default();
        genotype_log_likelihood_row(
            &StutterSubstitutionEmission,
            &evidence,
            SsrLocusParameters {
                candidates: &candidates,
                contexts: SsrScoringContextTable::new(&contexts, candidates.len()),
                outlier_weight: DEFAULT_OUTLIER_WEIGHT,
                reachable_length_count: 12,
            },
            &view,
            &mut out,
            &mut scratch,
        );
    }

    /// **A sample with no observations at this tract leaves every genotype at zero**, without a
    /// branch: the sum is empty, and the prior is what decides. Pinned because an empty row is
    /// the shape a caller is most likely to hand in without meaning to.
    #[test]
    fn a_sample_with_no_reads_leaves_every_genotype_at_zero() {
        let fixture = Fixture::of(b"CA", &[4, 5]);
        let row = score_row(&fixture, &[], 2, 12);
        assert!(
            row.iter().all(|slot| slot.0 == 0.0),
            "an empty row should be all zeros, not {row:?}"
        );
    }
}
