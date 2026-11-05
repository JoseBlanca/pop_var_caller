//! **Which other samples a sample copies from in a window** (`doc/devel/ng/spec/neighbour_evidence_pass.md`
//! §3.3).
//!
//! **Disagreement against what the cohort's frequencies predict.** At each record where the sample
//! and a candidate are both informative, the observed disagreement is the chance that one allele
//! drawn from each differs, `1 − Σ_a q_self(a) · q_other(a)`; the expected one is the same against
//! the cohort's frequencies, `1 − Σ_a q_self(a) · p(a)` — what an unrelated sample would show.
//! Observed `d` out of expected `e`, summed over some records, is a rate, scored pessimistically:
//! the mean of its Beta(`d + 1`, `e − d + 1`) posterior plus two of its standard deviations,
//! [`disagreement_bound`]. A record where the whole cohort carries one allele says nothing of who is
//! related and counts nothing (its `e` is 0); a candidate whose shares are uncertain — near the
//! frequencies — scores as an unrelated one, not a close one (owner, 2026-11-05).
//!
//! **Judged stretch by stretch, ranked by stretches won** (owner, 2026-11-05). The compared records
//! are cut, in order, into stretches that each hold [`STRETCH_INFORMATION`] of expected
//! disagreement — about two common SNPs. In each stretch the [`STRETCH_WINNERS`] candidates with the
//! lowest scores there win it, with any tied with the last of them; a candidate with nothing to be
//! judged on in a stretch cannot win it. Candidates rank by stretches won, then by their score over
//! the whole window, then by sample order, and the first `K` are the templates. A template's
//! stretches won also set how often a copy that changes template lands on it (§3.3, §3.5).
//!
//! Why stretches: a sample's two copies change template along the window at every crossover in its
//! ancestry. A whole-window score ranks highest whoever is fairly close everywhere; the samples a
//! copy actually copies are exactly close over a part each, and only stretches short enough to fall
//! between crossovers see that. On a simulated F2 at 1× the whole-window rule left an offspring's
//! parents out of its 16 templates or among siblings that blurred them; stretch wins raised the
//! heterozygous offspring genotypes called right from 70 to 81 in 100
//! (`doc/devel/reports/implementations/local_sources_2026-11-05.md`).
//!
//! **Data wins near ties, and is never penalised for being certain.** Within a stretch the
//! pessimistic score charges a candidate for the little it can be judged on, so a deep sample beats a
//! shallow one that agrees as well — and it is the deep sample whose own data fills a site the target
//! has none at. A shallow candidate is judged only where it has reads, so a record without them is
//! not a disagreement; but it cannot win a stretch it has no reads in, so a shallow close relative
//! wins fewer stretches than a deep one would, and a deep sample that is fairly close may rank above
//! it.
//!
//! The comparison runs over at most [`AGREEMENT_RECORDS`] of the sample's informative records,
//! spread evenly by record order, so the choice costs the same in a dense window as in a sparse one
//! and is the same on every run. The background template is not chosen here; it joins every walk
//! (§3.3).

use super::model::GenotypePairs;
use super::records::{SampleRow, SiteRecord};
use std::ops::Range;

/// The most templates a walk copies from, the background aside (`K`, §3.3) — what the prototype
/// ran; fewer measured worse on every cohort tried (2, 4 and 8, 2026-11-05).
pub const MAX_TEMPLATES: usize = 16;

/// The most of a sample's informative records its templates are chosen over (§3.3).
pub const AGREEMENT_RECORDS: usize = 500;

/// **How much a stretch must hold before its candidates are judged in it** (§3.3): a stretch of
/// compared records ends once the disagreement an unrelated sample would show the sample over it,
/// `Σ 1 − Σ_a q_self(a) · p(a)`, reaches this. One is about two SNPs at frequency ½ where the sample
/// is certain, or twenty at 5%. Larger stretches measured worse (2026-11-05).
pub const STRETCH_INFORMATION: f64 = 1.0;

/// **How many candidates win each stretch**, besides any tied with the last of them (§3.3): one a
/// copy. Four measured worse than two (2026-11-05).
pub const STRETCH_WINNERS: usize = 2;

/// **A template's landing weight before the stretches it won** (§3.3, §3.5): its weight is
/// stretches won plus this, so a template that won none is still landed on, half as often as one
/// that won a single stretch.
pub const LANDING_WEIGHT_BASE: f64 = 0.5;

/// **One chosen template**: who it is — a sample, or a sample's haplotype — and how many stretches
/// of the window it won (§3.3).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ChosenTemplate<T> {
    /// The sample, or the `(sample, copy)` haplotype.
    pub source: T,
    /// How many of the window's stretches it was among the winners of.
    pub stretches_won: usize,
}

impl<T> ChosenTemplate<T> {
    /// **How likely, relative to the other templates, a copy that changes template is to land on
    /// this one** (§3.3, §3.5): stretches won plus [`LANDING_WEIGHT_BASE`].
    #[must_use]
    pub fn landing_weight(&self) -> f64 {
        self.stretches_won as f64 + LANDING_WEIGHT_BASE
    }
}

/// **The templates sample `sample` copies from in the window `records`**, best first:
/// `min(max_templates, samples − 1)` other samples, each with the stretches it won. `None` when the
/// sample has no informative record in the window — then there is nothing to choose by, and the
/// sample is not re-called there (§3.3).
#[must_use]
pub fn choose_templates(
    records: &[SiteRecord],
    sample: usize,
    max_templates: usize,
) -> Option<Vec<ChosenTemplate<usize>>> {
    let compared = compared_records(records, sample);
    let &first = compared.first()?;
    let samples = records[first].sample_count();
    let others: Vec<usize> = (0..samples).filter(|&other| other != sample).collect();
    let expected = compared
        .iter()
        .map(|&index| {
            let record = &records[index];
            allele_mismatch(record.template_shares(sample), record.background_shares())
        })
        .collect();
    let observed = others
        .iter()
        .flat_map(|&other| {
            compared.iter().map(move |&index| {
                let record = &records[index];
                (record.row(other) == SampleRow::Informative).then(|| {
                    allele_mismatch(
                        record.template_shares(sample),
                        record.template_shares(other),
                    )
                })
            })
        })
        .collect();
    let ranked = Disagreements { expected, observed }.ranked();
    Some(
        ranked
            .into_iter()
            .take(max_templates)
            .map(|candidate| ChosenTemplate {
                source: others[candidate.index],
                stretches_won: candidate.stretches_won,
            })
            .collect(),
    )
}

/// **The records a sample's templates are chosen over**, as indices into `records`: its
/// informative sites, at most [`AGREEMENT_RECORDS`] of them spread evenly by record order (§3.3).
/// Empty when it has none.
#[must_use]
pub(crate) fn compared_records(records: &[SiteRecord], sample: usize) -> Vec<usize> {
    let informative: Vec<usize> = (0..records.len())
        .filter(|&index| {
            records[index].is_site() && records[index].row(sample) == SampleRow::Informative
        })
        .collect();
    if informative.len() <= AGREEMENT_RECORDS {
        informative
    } else {
        (0..AGREEMENT_RECORDS)
            .map(|draw| informative[draw * informative.len() / AGREEMENT_RECORDS])
            .collect()
    }
}

/// **How [`choose_haplotypes`] reads a haplotype**: `(other sample, copy, record, out)` writes that
/// copy's allele probabilities at the record into `out`, which it clears first.
pub(crate) type HaplotypeReader<'a> = dyn Fn(usize, usize, usize, &mut Vec<f64>) + 'a;

/// **The haplotypes sample `sample` copies from in a phasing round** (phasing spec §2.3), best
/// first, as `(sample, copy)` with the stretches each won: at most `max_templates`, at most one of
/// each other sample.
///
/// Judged as [`choose_templates`] judges samples, stretch by stretch, with the disagreement read
/// off the sample's genotype from the previous round rather than its shares: at each record the
/// observed disagreement is the chance the haplotype carries an allele that genotype lacks,
/// `Σ_a h(a) · P(genotype without a)`, and the expected one the same for an allele drawn from the
/// cohort's frequencies, `Σ_a p(a) · P(genotype without a)` — `h` the haplotype, the genotype
/// probabilities `previous_posteriors`, one per `compared` record in its order. Both haplotypes of
/// every other sample compete in each stretch; each sample offers the one that won more stretches
/// (on a tie, the better whole-window score; then copy 0), and the offers rank as samples do in
/// [`choose_templates`], then by copy. `haplotype(other, copy, record, out)` writes the previous
/// round's haplotype into `out`, which it clears first.
///
/// One of each sample because without it the slots fill with copies of the commonest haplotype:
/// 871 wrong genotypes against 689 on the prototype's first cohort (spec §2.3).
///
/// # Panics
///
/// If `previous_posteriors` is not one posterior a compared record.
#[must_use]
pub(crate) fn choose_haplotypes(
    records: &[SiteRecord],
    sample: usize,
    compared: &[usize],
    previous_posteriors: &[Vec<f64>],
    haplotype: &HaplotypeReader<'_>,
    max_templates: usize,
) -> Vec<ChosenTemplate<(usize, usize)>> {
    let Some(&first) = compared.first() else {
        return Vec::new();
    };
    assert_eq!(
        compared.len(),
        previous_posteriors.len(),
        "one previous posterior a compared record"
    );
    let samples = records[first].sample_count();
    let mut pairs = GenotypePairs::default();
    // P(genotype without a) at each compared record, for every allele.
    let lacking: Vec<Vec<f64>> = compared
        .iter()
        .zip(previous_posteriors)
        .map(|(&index, posterior)| {
            let alleles = records[index].allele_count;
            (0..alleles)
                .map(|allele| {
                    pairs
                        .of(alleles)
                        .iter()
                        .zip(posterior)
                        .filter(|&(&(a, b), _)| a != allele && b != allele)
                        .map(|(_, &probability)| probability)
                        .sum()
                })
                .collect()
        })
        .collect();
    let expected = compared
        .iter()
        .zip(&lacking)
        .map(|(&index, lacked)| {
            records[index]
                .background_shares()
                .zip(lacked)
                .map(|(frequency, lacked)| frequency * lacked)
                .sum()
        })
        .collect();
    // Candidate `2·j + copy` is copy `copy` of the `j`-th other sample.
    let others: Vec<usize> = (0..samples).filter(|&other| other != sample).collect();
    let mut observed = Vec::with_capacity(2 * others.len() * compared.len());
    let mut carried = Vec::new();
    for &other in &others {
        for copy in 0..2 {
            for (&index, lacked) in compared.iter().zip(&lacking) {
                observed.push(
                    (records[index].row(other) == SampleRow::Informative).then(|| {
                        haplotype(other, copy, index, &mut carried);
                        carried
                            .iter()
                            .zip(lacked)
                            .map(|(carried, lacked)| carried * lacked)
                            .sum()
                    }),
                );
            }
        }
    }
    let judged = Disagreements { expected, observed }.judged();
    let mut offers: Vec<Judged> = judged
        .as_chunks::<2>()
        .0
        .iter()
        .map(|copies| {
            if copies[1].ranks_before(&copies[0]) {
                copies[1]
            } else {
                copies[0]
            }
        })
        .collect();
    offers.sort_by(Judged::order);
    offers
        .into_iter()
        .take(max_templates)
        .map(|offer| ChosenTemplate {
            source: (others[offer.index / 2], offer.index % 2),
            stretches_won: offer.stretches_won,
        })
        .collect()
}

/// **Every candidate's disagreement with the sample at each compared record**, observed against
/// expected (§3.3).
struct Disagreements {
    /// What an unrelated sample would show at each compared record — the same for every candidate.
    expected: Vec<f64>,
    /// What each candidate shows, candidate by candidate, one entry a compared record; `None` where
    /// the candidate is not informative.
    observed: Vec<Option<f64>>,
}

/// **One candidate as the stretches judged it**: its index among the candidates, the stretches it
/// won, and its score over the whole window.
#[derive(Clone, Copy, Debug)]
struct Judged {
    index: usize,
    stretches_won: usize,
    score: f64,
}

impl Judged {
    /// More stretches won first, then the lower whole-window score, then the earlier candidate.
    fn order(one: &Self, other: &Self) -> std::cmp::Ordering {
        other
            .stretches_won
            .cmp(&one.stretches_won)
            .then(one.score.total_cmp(&other.score))
            .then(one.index.cmp(&other.index))
    }

    fn ranks_before(&self, other: &Self) -> bool {
        Self::order(self, other).is_lt()
    }
}

impl Disagreements {
    fn records(&self) -> usize {
        self.expected.len()
    }

    /// Candidate `candidate`'s observed disagreement at each record.
    fn of(&self, candidate: usize) -> &[Option<f64>] {
        let records = self.records();
        &self.observed[candidate * records..(candidate + 1) * records]
    }

    /// **A candidate's pessimistic score over `stretch`**: [`disagreement_bound`] of its observed
    /// over its expected disagreement, summed where it is informative, the observed capped at the
    /// expected — a candidate less alike than an unrelated sample scores as one. Infinite when
    /// nothing was expected there: no shared record, or only records every sample agrees at.
    fn score(&self, candidate: usize, stretch: Range<usize>) -> f64 {
        let (observed, expected) = self.of(candidate)[stretch.clone()]
            .iter()
            .zip(&self.expected[stretch])
            .filter_map(|(observed, &expected)| observed.map(|observed| (observed, expected)))
            .fold((0.0, 0.0), |(o, e), (observed, expected)| {
                (o + observed, e + expected)
            });
        if expected > 0.0 {
            disagreement_bound(f64::min(observed, expected), expected)
        } else {
            f64::INFINITY
        }
    }

    /// **The stretches the compared records are cut into**, in order: each ends once the expected
    /// disagreement over it reaches [`STRETCH_INFORMATION`], counting only records where some
    /// candidate is informative; the last holds what is left.
    fn stretches(&self) -> Vec<Range<usize>> {
        let candidates = self.observed.len() / self.records().max(1);
        let mut stretches = Vec::new();
        let (mut start, mut held) = (0, 0.0);
        for record in 0..self.records() {
            if (0..candidates).any(|candidate| self.of(candidate)[record].is_some()) {
                held += self.expected[record];
            }
            if held >= STRETCH_INFORMATION {
                stretches.push(start..record + 1);
                (start, held) = (record + 1, 0.0);
            }
        }
        if start < self.records() {
            stretches.push(start..self.records());
        }
        stretches
    }

    /// **Every candidate, with the stretches it won and its whole-window score**, in candidate
    /// order. A stretch's winners are the [`STRETCH_WINNERS`] lowest scores there and any tied with
    /// the last of them, among the candidates with a finite score.
    fn judged(&self) -> Vec<Judged> {
        let candidates = self.observed.len() / self.records().max(1);
        let mut stretches_won = vec![0_usize; candidates];
        let mut scores = Vec::with_capacity(candidates);
        for stretch in self.stretches() {
            scores.clear();
            scores.extend(
                (0..candidates)
                    .map(|candidate| self.score(candidate, stretch.clone()))
                    .enumerate()
                    .filter(|(_, score)| score.is_finite()),
            );
            if scores.is_empty() {
                continue;
            }
            let last = STRETCH_WINNERS.min(scores.len()) - 1;
            let (_, &mut (_, threshold), _) =
                scores.select_nth_unstable_by(last, |one, other| one.1.total_cmp(&other.1));
            for &(candidate, score) in &scores {
                if score <= threshold {
                    stretches_won[candidate] += 1;
                }
            }
        }
        stretches_won
            .into_iter()
            .enumerate()
            .map(|(index, stretches_won)| Judged {
                index,
                stretches_won,
                score: self.score(index, 0..self.records()),
            })
            .collect()
    }

    /// [`Self::judged`], best first ([`Judged::order`]).
    fn ranked(&self) -> Vec<Judged> {
        let mut judged = self.judged();
        judged.sort_by(Judged::order);
        judged
    }
}

/// **A pessimistic disagreement rate** from `disagreement` out of `shared`, the disagreement an
/// unrelated sample would have shown: the mean of the rate's Beta(`disagreement + 1`,
/// `shared − disagreement + 1`) posterior, from a flat prior, plus two of its standard deviations.
/// With nothing shared it is ½ + 2·√(1/12), above any rate measured.
#[must_use]
pub fn disagreement_bound(disagreement: f64, shared: f64) -> f64 {
    let mean = (disagreement + 1.0) / (shared + 2.0);
    mean + 2.0 * (mean * (1.0 - mean) / (shared + 3.0)).sqrt()
}

/// **The chance that one allele drawn from each of two allele-share vectors differs**:
/// `1 − Σ_a p(a) · q(a)`. Two certain and equal samples give 0; a sample whose shares are the
/// cohort's frequencies gives, against any other, what an unrelated sample would.
fn allele_mismatch(one: impl Iterator<Item = f64>, other: impl Iterator<Item = f64>) -> f64 {
    1.0 - one.zip(other).map(|(p, q)| p * q).sum::<f64>()
}

#[cfg(test)]
mod tests;
