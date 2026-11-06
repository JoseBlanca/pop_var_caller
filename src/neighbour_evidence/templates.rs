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
//! are cut, in order, into stretches that each hold [`stretch_information`] of expected
//! disagreement — one unit, about two common SNPs, up to 19 other samples, and more in a larger
//! cohort, so that an unrelated sample matching a stretch by chance does not take a relative's
//! place (§3.3). In each stretch the [`STRETCH_WINNERS`] candidates with the
//! lowest scores there win it, with any tied with the last of them; a candidate with nothing to be
//! judged on in a stretch cannot win it. Candidates rank by stretches won, then by their score over
//! the whole window, then by sample order, and the first `K` are the templates. A template's
//! stretches won also set how often a copy that changes template lands on it (§3.3, §3.5).
//!
//! Why stretches: a sample's two copies change template along the window at every crossover in its
//! ancestry. A whole-window score ranks highest whoever is fairly close everywhere; the samples a
//! copy actually copies are exactly close over a part each, and only stretches short enough to fall
//! between crossovers see that. On a simulated F2 at 1× the whole-window rule left an offspring's
//! parents out of its 16 templates or among siblings that blurred them. On two simulated F2s (two
//! inbred parents at 30×, 20 offspring at 1×), stretch wins raised the heterozygous offspring
//! genotypes called right from 72 and 67 in 100 to 82 and 81
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
/// ran; 2, 4 and 8 measured worse on every cohort of the 2026-11-05 report (40 samples at 1× to 30×,
/// and an F2 of 20 offspring at 1×).
pub const MAX_TEMPLATES: usize = 16;

/// The most of a sample's informative records its templates are chosen over (§3.3).
pub const AGREEMENT_RECORDS: usize = 500;

/// **The least a stretch holds before its candidates are judged in it** (§3.3), in units of the
/// disagreement an unrelated sample would show the sample over it, `Σ 1 − Σ_a q_self(a) · p(a)`.
/// One is about two SNPs at frequency ½ where the sample is certain, or twenty at 5%. Small cohorts
/// hold their stretches to this: on an F2 of 20 offspring at 1×, two units called 3 to 13 more
/// heterozygotes in 1,000 wrong or missing than one did, with crossovers as dense as the
/// benchmark's or ten times rarer (`stretch_scale_2026-11-06.md`).
pub const MIN_STRETCH_INFORMATION: f64 = 1.0;

/// **How many chance matches a stretch may draw** (§3.3), the scale of [`stretch_information`]: a
/// stretch holds `ln(candidates / this)` units, so that — if an unrelated candidate matches a
/// stretch of `I` units as well as a relative with chance about `e^-I` — about this many of the
/// cohort's candidates do. Measured, not derived: the model is crude, and 7 is the value that holds
/// an F2 of 20 offspring at one unit while a 400-sample outbred cohort gets four, which recovered
/// the phase that one unit lost there (98,871 → 108,600 heterozygotes phased right, of 145,951)
/// and called 570 fewer true variants wrong (`stretch_scale_2026-11-06.md`).
pub const TOLERATED_CHANCE_MATCHES: f64 = 7.0;

/// **How many samples win each stretch**, besides any tied with the last of them (§3.3): one a copy.
/// Four measured worse than two on the cohorts of the 2026-11-05 report.
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
    let mut sums = StretchSums::new(
        expected,
        &judged_records(records, sample, &compared),
        stretch_information(others.len()),
        1,
    );
    for &other in &others {
        sums.add(compared.iter().map(|&index| {
            let record = &records[index];
            (record.row(other) == SampleRow::Informative).then(|| {
                allele_mismatch(
                    record.template_shares(sample),
                    record.template_shares(other),
                )
            })
        }));
    }
    let ranked = sums.ranked();
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
/// every other sample compete in each stretch, a sample taking at most one of its two winners'
/// places (so a homozygous sample's two equal copies do not take both); each sample offers the one
/// that won more stretches
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
    let mut sums = StretchSums::new(
        expected,
        &judged_records(records, sample, compared),
        stretch_information(others.len()),
        2,
    );
    let mut carried = Vec::new();
    for &other in &others {
        for copy in 0..2 {
            sums.add(compared.iter().zip(&lacking).map(|(&index, lacked)| {
                (records[index].row(other) == SampleRow::Informative).then(|| {
                    haplotype(other, copy, index, &mut carried);
                    carried
                        .iter()
                        .zip(lacked)
                        .map(|(carried, lacked)| carried * lacked)
                        .sum()
                })
            }));
        }
    }
    let judged = sums.judged();
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

/// **Whether some candidate is judged at each compared record**: whether any sample but `sample`
/// is informative there — what a stretch's information counts (§3.3).
fn judged_records(records: &[SiteRecord], sample: usize, compared: &[usize]) -> Vec<bool> {
    compared
        .iter()
        .map(|&index| {
            let record = &records[index];
            (0..record.sample_count())
                .any(|other| other != sample && record.row(other) == SampleRow::Informative)
        })
        .collect()
}

/// **How much a stretch holds when `candidates` other samples compete for it** (§3.3):
/// `ln(candidates / TOLERATED_CHANCE_MATCHES)` units, never less than
/// [`MIN_STRETCH_INFORMATION`] — one unit up to 19 candidates, 1.7 at 39, 4.0 at 399, 4.9 at 999
/// and 6.6 at 4,999. Short stretches see a copy's change of source where it happens, but among
/// many candidates some unrelated one matches a short stretch by chance; at 400 samples one unit
/// gave a stretch 6.5 winners where two were asked for.
#[must_use]
pub fn stretch_information(candidates: usize) -> f64 {
    (candidates as f64 / TOLERATED_CHANCE_MATCHES)
        .ln()
        .max(MIN_STRETCH_INFORMATION)
}

/// **The stretches the compared records are cut into**, in order (§3.3): each ends at the record
/// where the expected disagreement summed over it reaches `information`, counting only records
/// where `judged[record]`. What is left after the last such record joins the stretch before it, so
/// no stretch is judged on less; a window holding less than `information` is one stretch.
fn stretches(expected: &[f64], judged: &[bool], information: f64) -> Vec<Range<usize>> {
    let mut stretches: Vec<Range<usize>> = Vec::new();
    let (mut start, mut held) = (0, 0.0);
    for (record, (&expected, &judged)) in expected.iter().zip(judged).enumerate() {
        if judged {
            held += expected;
        }
        if held >= information {
            stretches.push(start..record + 1);
            (start, held) = (record + 1, 0.0);
        }
    }
    match stretches.last_mut() {
        Some(last) => last.end = expected.len(),
        None if !expected.is_empty() => stretches.push(0..expected.len()),
        None => {}
    }
    stretches
}

/// **Every candidate's disagreement with the sample, summed stretch by stretch** (§3.3): observed
/// against expected, summed over the records of each stretch where the candidate is informative.
/// Only the sums are kept — a pair a candidate a stretch, not a number a candidate a record.
struct StretchSums {
    /// Each compared record's stretch.
    stretch_of: Vec<usize>,
    /// How many stretches there are.
    stretches: usize,
    /// What an unrelated sample would show at each compared record — the same for every candidate.
    expected: Vec<f64>,
    /// How many candidates are one sample's: 1 when candidates are samples, 2 when they are a
    /// sample's two haplotypes, which are consecutive.
    copies: usize,
    /// Each candidate's `(observed, expected)` in each stretch, candidate by candidate.
    sums: Vec<(f64, f64)>,
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

impl StretchSums {
    /// Stretches of `information` cut over `expected` where `judged`, no candidate added yet;
    /// `copies` candidates to a sample.
    fn new(expected: Vec<f64>, judged: &[bool], information: f64, copies: usize) -> Self {
        let cut = stretches(&expected, judged, information);
        let mut stretch_of = vec![0; expected.len()];
        for (stretch, range) in cut.iter().enumerate() {
            stretch_of[range.clone()].fill(stretch);
        }
        Self {
            stretch_of,
            stretches: cut.len(),
            expected,
            copies,
            sums: Vec::new(),
        }
    }

    /// Adds the next candidate: its observed disagreement at each compared record, `None` where it
    /// is not informative.
    fn add(&mut self, observed: impl IntoIterator<Item = Option<f64>>) {
        let start = self.sums.len();
        self.sums.resize(start + self.stretches, (0.0, 0.0));
        let sums = &mut self.sums[start..];
        for ((observed, &stretch), &expected) in observed
            .into_iter()
            .zip(&self.stretch_of)
            .zip(&self.expected)
        {
            if let Some(observed) = observed {
                sums[stretch].0 += observed;
                sums[stretch].1 += expected;
            }
        }
    }

    fn candidates(&self) -> usize {
        self.sums.len() / self.stretches.max(1)
    }

    /// **Every candidate, with the stretches it won and its whole-window score**, in candidate
    /// order. A stretch's winners are the candidates of the [`STRETCH_WINNERS`] samples with the
    /// lowest scores there — a sample scoring by its better copy — and any tied with the last of
    /// them, among the candidates with a finite score: so a sample whose two haplotypes are the
    /// same cannot take both places.
    fn judged(&self) -> Vec<Judged> {
        let candidates = self.candidates();
        let samples = candidates / self.copies;
        let mut stretches_won = vec![0_usize; candidates];
        let mut scores = vec![f64::INFINITY; candidates];
        let mut best = Vec::with_capacity(samples);
        for stretch in 0..self.stretches {
            for (candidate, score) in scores.iter_mut().enumerate() {
                *score = pessimistic_score(self.sums[candidate * self.stretches + stretch]);
            }
            best.clear();
            best.extend(
                scores
                    .chunks_exact(self.copies)
                    .map(|copies| copies.iter().copied().fold(f64::INFINITY, f64::min))
                    .filter(|score| score.is_finite()),
            );
            if best.is_empty() {
                continue;
            }
            let last = STRETCH_WINNERS.min(best.len()) - 1;
            let (_, &mut threshold, _) = best.select_nth_unstable_by(last, f64::total_cmp);
            for (won, &score) in stretches_won.iter_mut().zip(&scores) {
                if score <= threshold {
                    *won += 1;
                }
            }
        }
        stretches_won
            .into_iter()
            .enumerate()
            .map(|(index, stretches_won)| {
                let whole = self.sums[index * self.stretches..(index + 1) * self.stretches]
                    .iter()
                    .fold((0.0, 0.0), |(o, e), &(observed, expected)| {
                        (o + observed, e + expected)
                    });
                Judged {
                    index,
                    stretches_won,
                    score: pessimistic_score(whole),
                }
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

/// **A candidate's pessimistic score from its summed `(observed, expected)`**:
/// [`disagreement_bound`], the observed capped at the expected — a candidate less alike than an
/// unrelated sample scores as one. Infinite when nothing was expected: no shared record, or only
/// records every sample agrees at.
fn pessimistic_score((observed, expected): (f64, f64)) -> f64 {
    if expected > 0.0 {
        disagreement_bound(f64::min(observed, expected), expected)
    } else {
        f64::INFINITY
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
