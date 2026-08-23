//! The comparator: Hardy–Weinberg at one estimated allele frequency, plugged in as though it
//! were the truth.
//!
//! It is the route this caller does *not* take, kept behind the same seam only so the change the
//! default makes stays measurable — **never a shipping default**
//! (`doc/devel/ng/arch/calling_priors.md` §7).
//!
//! ## What it does differently, and it is one line of algebra
//!
//! Both implementations behind the seam turn a concentration — chromosomes the prior behaves as
//! though it had already seen — into one log-probability per candidate genotype, and both apply
//! the same two-branch inbreeding mixture (`doc/devel/ng/spec/calling_priors.md` §3.2). They
//! differ only in the random-mating branch:
//!
//! - the **default** averages the genotype probability over every allele frequency the
//!   concentration finds plausible, weighting each by how plausible it is
//!   ([`dirichlet_multinomial`](super::dirichlet_multinomial));
//! - **this one** collapses the concentration to a single frequency, `α_a / Σα`, and evaluates
//!   the genotype probability there.
//!
//! Genotype probability is quadratic in the frequency, and the average of a curve is not the
//! curve at the average. For a homozygote `E[p²] = p̄² + Var(p)`, so **plugging in undercounts
//! homozygotes by exactly the variance of the frequency** — the same term for both homozygotes —
//! and since the row sums to one that mass lands on the heterozygotes, which receive twice what
//! either homozygote loses (spec §2.2).
//!
//! **`Var(p)` is how badly the frequency is pinned down.** With a thousand samples it is
//! negligible and the two implementations agree; with one sample at low depth it dominates. So
//! the gap is largest precisely in the corner this caller commits to supporting, which is why the
//! comparator exists rather than the difference being argued about.
//!
//! ## The measurement it is kept for
//!
//! On the GIAB trio, each sample called on its own at 5×, swapping this route for the default
//! took SNP genotype accuracy at true variants from **83.6% to 94.6%**, and the sites where a
//! sample carrying two copies of the variant was called heterozygous from **214 to 8**, with the
//! emitted variant set byte-identical — the prior moves genotypes among emitted variants and does
//! not change which sites are emitted (spec §2.2). That is one corner: one sample at a time,
//! high-quality human data, 5×, and it is the corner the change was aimed at.
//!
//! ## The trap this file must not fall into
//!
//! **The gain in that measurement is not "marginalize". It is the starting concentration**
//! (spec §2.3). Production's plug-in path regularised its frequency estimate with a reference
//! pseudocount of 10, which put the prior odds on a heterozygote at 22:1 at the configuration
//! that mattered; marginalizing over that same concentration gives 20:1 — the same wrong answer,
//! computed more expensively. **So the comparator runs on the same concentration as the default
//! and supplies no pseudocount of its own**, and
//! [`tests::the_row_is_hardy_weinberg_at_the_handed_concentration_and_nothing_else`] is what holds
//! that: it checks the row against a closed form evaluated at exactly the frequencies handed in,
//! which any hidden pseudocount would break.

use crate::genetics::PROBABILITY_FLOOR;
use crate::ng::calling::genotype_prior::dirichlet_multinomial::log_sum_exp_2;
use crate::ng::calling::genotype_prior::{GenotypePriorModel, PriorRow};
use crate::ng::types::{InbreedingF, LogProb};

/// The comparator implementation of the step-8 seam: Hardy–Weinberg at the plug-in frequency
/// `α'_s(a) / Σα'_s`, with the same inbreeding mixture as the default.
///
/// **Kept only for the spec's change measurements and the production differential.** A run that
/// selects it is measuring the prior, not calling variants.
///
/// **The derives match its sibling's**, for its sibling's reason: a run that compares two priors
/// has to be able to say which arm produced a row.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct PlugInWrightPrior;

impl GenotypePriorModel for PlugInWrightPrior {
    fn name(&self) -> &'static str {
        "plug-in-wright"
    }

    fn fill_genotype_log_priors(&self, row: &mut PriorRow<'_>, inbreeding: InbreedingF) {
        fill_plug_in_mixture_log_priors(row, inbreeding.get());
    }
}

/// The mixture on a bare coefficient rather than on the newtype, for the reason its twin in
/// [`dirichlet_multinomial`](super::dirichlet_multinomial) has one: `F = 1` is the mathematical
/// edge of the model, and each implementation is pinned there **against its own closed form**.
///
/// **Not against each other**, and a test that tried would fail for a legitimate reason: at
/// `F = 1` the two still differ by `ln((Σα + 1)/Σα)` — 0.69 nats at a single sample's seed —
/// because complete inbreeding scales the random-mating branch rather than removing it, so the
/// genuine gap between averaging `pq` and squaring its average survives.
///
/// ```text
/// p_a           = α_a / Σα                                   the plug-in frequency
/// log P_random  = ln C(ploidy; counts) + Σ_a counts_a · ln p_a
/// homozygous g  = logsumexp( ln(1 − F) + log P_random,  ln F + ln p_a )
/// any other g   = ln(1 − F) + log P_random
/// ```
///
/// **The random-mating branch is already a true probability here**, and that is the one
/// structural difference from the default's mixture. There, the ported primitive drops the
/// genotype-independent `lgamma(Σα + m) − lgamma(Σα)` because it cancels within a row, and the
/// mixture has to put it back before mixing an unnormalised branch against a normalised one. A
/// multinomial at fixed frequencies drops nothing, so there is no term to restore — **and a
/// reader who ports this back the other way must not add one**.
///
/// **At two alleles and two copies, and at any `F` below 1, this is the textbook Wright
/// formula** — which is what makes `wright_genotype_log_priors` (`src/genetics.rs`) an oracle for
/// it rather than a second code path: `(1 − F)p² + Fp` is `p² + Fpq`, and `(1 − F)·2pq` is
/// `2pq(1 − F)`. The mixture form is the one built because it keeps saying something at three
/// alleles and at a tetraploid, which the biallelic formulas do not (spec §3.2).
///
/// **The two part company at `F = 1` exactly, and only there, because they floor different
/// things.** Production floors the finished probability, so a heterozygote whose probability is
/// exactly zero lands at `ln(1e-300) = −690.78`; this floors the *weight* `1 − F` and then adds
/// the genotype's own `ln 2pq`, landing 6.73 nats lower at `−697.50`. Both say "impossible" and
/// neither says `−∞`, which is the contract; they disagree about how impossible. At `F = 0.99`
/// they still agree to 9e-16, so the divergence is the endpoint and nothing before it — which is
/// why the oracle test sweeps to `0.9` and the `F = 1` limit is pinned separately, against the
/// closed form rather than against production.
///
/// **`1 − F` is floored and `F` is not**, matching the default exactly, and the two now share
/// [`log_sum_exp_2`] so they cannot drift apart about `−∞`. The frequencies are floored for the
/// same reason: a genotype the prior rules out carries a finite, very negative log-prior rather
/// than `−∞` (spec §8), which is also what production's own biallelic version does internally.
///
/// ## Cost
///
/// One `ln` per allele and one multiply-add per (genotype, non-zero allele) pair, plus one
/// `logsumexp` per homozygous genotype. **No `lgamma` at all**, which is the whole of the
/// default's expense — so the comparator is cheaper than the thing it is compared against, and a
/// measurement that finds it faster has found nothing.
fn fill_plug_in_mixture_log_priors(row: &mut PriorRow<'_>, inbreeding: f64) {
    debug_assert!(
        (0.0..=1.0).contains(&inbreeding),
        "the inbreeding coefficient must be a fraction in [0, 1]; got {inbreeding}. A value \
         outside it — or a NaN — makes the whole row NaN rather than failing, and a negative one \
         poisons only the homozygotes, which normalises to a plausible wrong answer"
    );

    let concentration = row.concentration().get();
    let allele_count = row.concentration().allele_count();
    let genotype_allele_counts = row.genotype_allele_counts();
    let log_multinomial_coeffs = row.log_multinomial_coeffs();
    let homozygous_allele_for = row.homozygous_allele_for();
    let concentration_total: f64 = concentration.iter().sum();

    // Floored for the same reason the default floors `1 − F`, and to the same constant: a
    // genotype the prior rules out gets a finite, very negative log-prior rather than `−∞`.
    let log_weight_identical_by_descent = inbreeding.ln();
    let log_weight_independent_draws = (1.0 - inbreeding).max(PROBABILITY_FLOOR).ln();

    let (log_frequency, out) = row.scratch_and_out();
    for (slot, &alpha) in log_frequency.iter_mut().zip(concentration) {
        *slot = (alpha / concentration_total).max(PROBABILITY_FLOOR).ln();
    }

    // `zip` would truncate to the shorter of the two, which is the silent failure this module's
    // own documentation names as the worst available here — `PriorRow::new` holds their equality
    // in release so it cannot happen.
    for (genotype, (slot, homozygous)) in out.iter_mut().zip(homozygous_allele_for).enumerate() {
        let counts = &genotype_allele_counts[genotype * allele_count..][..allele_count];
        let mut log_random_mating = log_multinomial_coeffs[genotype];
        for (&copies, &log_p) in counts.iter().zip(log_frequency.iter()) {
            // Skipping the zero counts is the same saving the default's primitive makes, and it
            // is a saving only: the frequencies are floored, so `0 × ln p` would be an ordinary
            // zero rather than a `NaN`.
            if copies > 0 {
                log_random_mating += f64::from(copies) * log_p;
            }
        }
        let independent_draws = log_weight_independent_draws + log_random_mating;
        *slot = LogProb(match homozygous {
            // The homozygous test is this lookup and nothing else — no comparison of the copy
            // counts anywhere in this function. That is what gives the above-diploidy question
            // one place to change (spec §3.3), and it is the same lookup the default consumes.
            Some(allele) => log_sum_exp_2(
                independent_draws,
                log_weight_identical_by_descent + log_frequency[usize::from(allele.0)],
            ),
            None => independent_draws,
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genetics::wright_genotype_log_priors;
    use crate::ng::calling::GenotypeTable;
    use crate::ng::calling::genotype_prior::{Concentration, MarginalizedDirichletPrior};
    use crate::ng::types::Ploidy;

    /// tomato1's fitted expected heterozygosity, spec §4.1 — the diversity every thin-input
    /// fixture here runs at, so a number quoted in one test means the same thing in another.
    const TOMATO_DIVERSITY: f64 = 6e-4;

    /// Run one implementation of the seam over one shape and hand back the row.
    ///
    /// The three genotypes of a biallelic diploid come back in the table's own order, which
    /// [`the_table_orders_the_biallelic_diploid_genotypes_as_the_oracle_does`] pins against the
    /// oracle's `(hom-ref, het, hom-alt)` rather than assuming.
    fn row_under(
        model: &impl GenotypePriorModel,
        concentration: &[f64],
        copies: u8,
        inbreeding: f64,
    ) -> Vec<f64> {
        let table = GenotypeTable::build(Ploidy::try_new(copies).unwrap(), concentration.len());
        let view = table.view();
        let mut scratch = vec![0.0; concentration.len()];
        let mut out = vec![LogProb(f64::NAN); view.genotype_count()];
        let mut row = PriorRow::new(
            Concentration::new(concentration),
            view.genotype_allele_counts(),
            view.log_multinomial_coeffs(),
            view.homozygous_alleles(),
            &mut scratch,
            &mut out,
        );
        model.fill_genotype_log_priors(&mut row, InbreedingF::try_new(inbreeding).unwrap());
        // **Checked here so no test below has to remember.** Every comparison in this file folds
        // its departures with `f64::max`, which *ignores* `NaN` — so a row whose entries were
        // never written would score a departure of zero and pass. Measured: filling only the
        // first genotype left `the_row_is_hardy_weinberg_at_the_handed_concentration_and_nothing_else`
        // green, and it is the test this module's own documentation names as its pin. The trait's
        // contract requires every entry finite, so asserting it here is free.
        assert!(
            out.iter().all(|entry| entry.get().is_finite()),
            "the seam's contract requires every entry finite, got {out:?}"
        );
        out.iter().map(|entry| entry.get()).collect()
    }

    /// A row shifted so its largest entry is zero — the comparison two implementations of this
    /// seam admit, because the contract is "up to a shared additive constant".
    fn shifted_to_peak(row: &[f64]) -> Vec<f64> {
        let peak = row.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        row.iter().map(|entry| entry - peak).collect()
    }

    /// A row turned into a genuine probability distribution — every entry divided by the total.
    ///
    /// **Both implementations need this before they can be compared entry by entry**, and only
    /// one of them needs it for a reason worth knowing: the primitive the default ports drops the
    /// genotype-independent term `lgamma(Σα + m) − lgamma(Σα)`, and the mixture restores it on
    /// **both** branches rather than on neither — so its row is a true distribution scaled by
    /// that constant. At the diploid `m = 2` every fixture here uses, the constant is
    /// `Σα(Σα + 1)`: measured, the default's raw row sums to 2.00180036 at `α = (1, 6e-4)`, and
    /// `1.0006 × 2.0006` is 2.00180036. The comparator's row is already a distribution.
    fn normalised(row: &[f64]) -> Vec<f64> {
        let total: f64 = row.iter().map(|entry| entry.exp()).sum();
        row.iter().map(|entry| entry.exp() / total).collect()
    }

    /// The genotype table's order for a biallelic diploid is `0/0`, `0/1`, `1/1`, which is the
    /// order [`wright_genotype_log_priors`] returns. **Pinned rather than assumed**, because every
    /// oracle comparison below indexes into the row by position.
    #[test]
    fn the_table_orders_the_biallelic_diploid_genotypes_as_the_oracle_does() {
        let table = GenotypeTable::build(Ploidy::try_new(2).unwrap(), 2);
        let view = table.view();
        assert_eq!(view.genotype_count(), 3);
        assert_eq!(view.genotype_allele_counts(), [2, 0, 1, 1, 0, 2]);
    }

    /// **The comparator is the textbook Wright formula wherever the textbook formula applies** —
    /// biallelic, diploid — which is what makes production's own
    /// [`wright_genotype_log_priors`] an independent oracle for it rather than a second code path.
    ///
    /// Swept over frequencies from 1 in 10,000 to a half and inbreeding from outbred to 0.9. The
    /// concentration is `(1 − p, p)` scaled by any total, because the plug-in route reads only the
    /// ratio — which this checks by running each case at two totals four orders apart and
    /// requiring the same row.
    #[test]
    fn the_comparator_is_the_wright_formula_at_two_alleles_and_two_copies() {
        let mut worst = 0.0_f64;
        for frequency in [1e-4, TOMATO_DIVERSITY, 1e-2, 0.1, 0.5] {
            for inbreeding in [0.0, 0.3, 0.8, 0.9] {
                let (hom_ref, het, hom_alt) = wright_genotype_log_priors(frequency, inbreeding);
                let oracle = shifted_to_peak(&[hom_ref, het, hom_alt]);
                for total in [1.0, 1e4] {
                    let concentration = [(1.0 - frequency) * total, frequency * total];
                    let row = shifted_to_peak(&row_under(
                        &PlugInWrightPrior,
                        &concentration,
                        2,
                        inbreeding,
                    ));
                    for (ours, theirs) in row.iter().zip(&oracle) {
                        worst = worst.max((ours - theirs).abs());
                    }
                }
            }
        }
        assert!(
            worst < 1e-12,
            "worst disagreement with the oracle: {worst:e}"
        );
    }

    /// **The row is Hardy–Weinberg at the concentration it was handed and at nothing else** — the
    /// pin spec §2.3 asks for, and the reason this file exists rather than a port of production's
    /// plug-in path.
    ///
    /// Production regularised its frequency estimate with a reference pseudocount of 10, which is
    /// where the defect the default repairs actually lived. A hidden pseudocount here would be
    /// invisible in a ratio test and invisible in the oracle test above, which supplies its own
    /// concentration — so this one builds the closed form from the **entries of the buffer that
    /// was passed in**, over a triallelic tetraploid where a wrong reference entry cannot hide.
    #[test]
    fn the_row_is_hardy_weinberg_at_the_handed_concentration_and_nothing_else() {
        let concentration = [1.0, TOMATO_DIVERSITY, 3.0 * TOMATO_DIVERSITY];
        let inbreeding = 0.3;
        let total: f64 = concentration.iter().sum();
        let frequency: Vec<f64> = concentration.iter().map(|alpha| alpha / total).collect();

        let table = GenotypeTable::build(Ploidy::try_new(4).unwrap(), concentration.len());
        let view = table.view();
        let row = row_under(&PlugInWrightPrior, &concentration, 4, inbreeding);

        let mut expected = Vec::with_capacity(view.genotype_count());
        for genotype in 0..view.genotype_count() {
            let counts = &view.genotype_allele_counts()[genotype * 3..][..3];
            let mut log_random = view.log_multinomial_coeffs()[genotype];
            for (&copies, &p) in counts.iter().zip(&frequency) {
                log_random += f64::from(copies) * p.ln();
            }
            let independent = (1.0 - inbreeding).ln() + log_random;
            expected.push(match view.homozygous_alleles()[genotype] {
                Some(allele) => {
                    let ibd = inbreeding.ln() + frequency[usize::from(allele.0)].ln();
                    let peak = independent.max(ibd);
                    peak + ((independent - peak).exp() + (ibd - peak).exp()).ln()
                }
                None => independent,
            });
        }

        // **Every frequency here is a fraction of the buffer's own total**, which is what makes a
        // hidden pseudocount visible. Measured by injecting production's reference pseudocount of
        // 10 into this fixture: the genotypes carrying **no** reference copy move by **9.2 nats**,
        // because every one of their frequencies is divided by a total ten times larger. That is
        // nine trillion times this tolerance.
        //
        // **The fully-reference genotype is the one that barely moves** — 0.0067 nats — since
        // `p_ref` goes from `1/1.0024` to `10/10.0024` and both are within a quarter of a per cent
        // of one. A reader checking this by eye should look at the genotypes without a reference
        // copy, not at the obvious one.
        let worst = row
            .iter()
            .zip(&expected)
            .map(|(ours, theirs)| (ours - theirs).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            worst < 1e-12,
            "worst departure from the closed form: {worst:e}"
        );
    }

    /// **The two implementations agree once the frequency is pinned down**, which is the whole
    /// content of the difference between them: they differ by the variance of the frequency, and
    /// that variance is `p(1 − p)/(Σα + 1)`.
    ///
    /// So scaling the concentration at a fixed ratio drives them together. Measured on a biallelic
    /// diploid at tomato's diversity, as the worst disagreement anywhere in the row once both are
    /// shifted to their own peak — the comparison the seam's contract admits, since a row is
    /// defined only up to a shared additive constant:
    ///
    /// ```text
    ///   Σα ×1      6.73 nats
    ///   Σα ×10²    2.86
    ///   Σα ×10⁴    0.154
    ///   Σα ×10⁶    0.00166
    /// ```
    ///
    /// **The last step closes the gap 93-fold and the first only 2.4-fold**, and that is the
    /// algebra rather than noise: the gap is dominated by the homozygote for the alternative
    /// allele, `ln(p̄²)` against `ln(p̄² + Var)`, which is only proportional to `Var` once
    /// `Var ≪ p̄²` — that is, once `Σα` is well past `1/p̄`, about 1,700 here. Below that the two
    /// priors are not near each other at all, which is the point of keeping the comparator.
    #[test]
    fn the_two_implementations_converge_as_the_frequency_becomes_certain() {
        let mut worst_at_scale = Vec::new();
        for scale in [1.0, 1e2, 1e4, 1e6] {
            let concentration = [scale, TOMATO_DIVERSITY * scale];
            let plug_in = shifted_to_peak(&row_under(&PlugInWrightPrior, &concentration, 2, 0.0));
            let marginalized = shifted_to_peak(&row_under(
                &MarginalizedDirichletPrior,
                &concentration,
                2,
                0.0,
            ));
            worst_at_scale.push(
                plug_in
                    .iter()
                    .zip(&marginalized)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0_f64, f64::max),
            );
        }

        for pair in worst_at_scale.windows(2) {
            assert!(pair[1] < pair[0], "{:?}", worst_at_scale);
        }
        // Once the variance is small against the squared frequency, a hundredfold more conviction
        // closes the gap about a hundredfold — the `1/(Σα + 1)` the algebra predicts.
        let closed_by = worst_at_scale[2] / worst_at_scale[3];
        assert!((90.0..100.0).contains(&closed_by), "closed by {closed_by}");
        assert!(worst_at_scale[3] < 2e-3, "{:e}", worst_at_scale[3]);
    }

    /// **At one sample's thin input the comparator undercounts homozygotes by exactly the
    /// variance of the frequency**, which is spec §2.2's algebra made a test rather than restated.
    ///
    /// The concentration is the seed a single sample starts from — the cohort term is exactly zero
    /// there (spec §6) — so the frequency is as badly pinned down as it ever gets. Both rows are
    /// normalised to sum to one first, because a row is defined only up to a shared additive
    /// constant and the default's is offset by `lgamma(Σα + m) − lgamma(Σα)`: mixing an
    /// unnormalised
    /// random-mating branch against a normalised inbreeding branch is the defect spec §3.2 records
    /// against production, and the default avoids it by offsetting **both** branches alike.
    ///
    /// **The identity is exact rather than asymptotic** — `E[p²] − p̄² = Var(p)` and
    /// `E[pq] − p̄q̄ = −Var(p)` are algebra — so the tolerance is 1e-9 and the worst measured
    /// residual is 5.3e-15. Measured, at tomato's diversity, `Var(p) = 2.9955e-4`:
    ///
    /// ```text
    ///                  default     comparator   difference
    ///   hom-ref        0.99910     0.99880      −Var
    ///   het            5.9946e-4   1.1986e-3    +2·Var
    ///   hom-alt        2.9991e-4   3.5957e-7    −Var
    /// ```
    ///
    /// So the comparator is **834 times less willing** to call a sample homozygous for the
    /// alternative allele, and hands that mass to the heterozygote — which it then calls at
    /// **3,333:1** against the homozygote where the default says **1.9988:1** — spec §2.3's
    /// `2·α_ref : (1 + α_alt)`, which is the 2:1 that section identifies as the defensible value,
    /// reached exactly only as the alternative concentration goes to zero.
    #[test]
    fn at_one_samples_thin_input_the_comparator_undercounts_homozygotes_by_the_variance() {
        let concentration = [1.0, TOMATO_DIVERSITY];
        let total: f64 = concentration.iter().sum();
        let frequency = concentration[1] / total;
        // A Dirichlet frequency's variance, `p(1 − p)/(Σα + 1)` — spec §1.
        let variance = frequency * (1.0 - frequency) / (total + 1.0);

        let plug_in = normalised(&row_under(&PlugInWrightPrior, &concentration, 2, 0.0));
        let default = normalised(&row_under(
            &MarginalizedDirichletPrior,
            &concentration,
            2,
            0.0,
        ));

        // Genotype 0 is hom-ref, 1 the heterozygote, 2 the homozygote for the alternative.
        for homozygote in [0, 2] {
            let taken = default[homozygote] - plug_in[homozygote];
            assert!(
                (taken / variance - 1.0).abs() < 1e-9,
                "homozygote {homozygote} lost {taken:e} against a variance of {variance:e}"
            );
        }
        let gained = plug_in[1] - default[1];
        assert!(
            (gained / (2.0 * variance) - 1.0).abs() < 1e-9,
            "the heterozygote gained {gained:e} against twice a variance of {variance:e}"
        );

        // The two sizes the spec quotes, against their closed forms rather than against a band —
        // including this one, which was a band against a hard-coded 833 and was the one number in
        // the block that turned out wrong. `E[p²]/p̄² = 1 + Var/p̄²`, which is 834.08 here.
        let undercount = 1.0 + variance / (frequency * frequency);
        assert!((default[2] / plug_in[2] / undercount - 1.0).abs() < 1e-9);
        assert!((834.0..834.2).contains(&undercount), "{undercount}");
        // Spec §2.3's ratio exactly: `2·α_ref : (1 + α_alt)`, which reaches 2:1 only as the
        // alternative concentration goes to zero and is 1.9988 here.
        let default_odds = 2.0 * concentration[0] / (1.0 + concentration[1]);
        assert!((default[1] / default[2] / default_odds - 1.0).abs() < 1e-9);
        assert!((1.9987..1.9989).contains(&default_odds), "{default_odds}");
        // Hardy-Weinberg's own `2q/p` for the comparator, which is 3,333 here.
        let plug_in_odds = 2.0 * (1.0 - frequency) / frequency;
        assert!((plug_in[1] / plug_in[2] / plug_in_odds - 1.0).abs() < 1e-9);
        assert!((3_333.0..3_334.0).contains(&plug_in_odds), "{plug_in_odds}");
    }

    /// The inbreeding mixture is the same one, so the comparator honours it exactly where the
    /// default does: at complete inbreeding every heterozygote sits at the probability floor and
    /// the two homozygotes stand in the ratio of their own frequencies.
    ///
    /// `InbreedingF::try_new(1.0)` returns `Ok` **today**; the prerequisites plan tightens the
    /// type to `[0, 1)`, and when it lands this test needs revisiting rather than deleting — the
    /// limit is worth pinning either way, through a raw-value path if the newtype stops admitting
    /// it (spec §7, §12 test 3).
    #[test]
    fn at_complete_inbreeding_the_comparator_leaves_only_the_homozygotes() {
        let concentration = [1.0, TOMATO_DIVERSITY];
        let total: f64 = concentration.iter().sum();
        let row = row_under(&PlugInWrightPrior, &concentration, 2, 1.0);

        // **The heterozygote is impossible but finite**, which is the contract: its only branch
        // is the one weighted `1 − F`, and that weight is floored rather than left at `ln(0)`. So
        // the entry is the floor's logarithm plus the genotype's own random-mating term — 697.5
        // below zero here, and more than 600 nats below either homozygote, which is impossible for
        // any purpose while still being a number.
        assert!(row[1].is_finite(), "{}", row[1]);
        assert!(
            row[1] < row[0] - 600.0 && row[1] < row[2] - 600.0,
            "{row:?}"
        );
        // The two homozygotes stand in the ratio of the frequencies themselves.
        let odds = (row[0] - row[2]).exp();
        assert!(
            (odds / (concentration[0] / concentration[1]) - 1.0).abs() < 1e-9,
            "{odds} against {}",
            concentration[0] / concentration[1]
        );
        // And with the heterozygote gone the two homozygotes carry the whole row, which is the
        // clearest statement that the mixture's inbreeding branch is a true probability.
        assert!((row.iter().map(|entry| entry.exp()).sum::<f64>() - 1.0).abs() < 1e-9);
        let _ = total;
    }

    /// **The floor on the frequencies is the only release-mode defence on this path**, and this
    /// is what pins it.
    ///
    /// Every value check in this folder is a `debug_assert!`, and `[profile.release]` sets
    /// `debug-assertions = false` — so in a shipping build a concentration entry of zero reaches
    /// this function unchecked, and the floor is the one thing standing between it and the `−∞`
    /// the seam's contract forbids by name.
    ///
    /// **It is defensive rather than live, and the size says so.** Deleting the floor changes no
    /// output until the ratio between an entry and the total underflows, which needs a total near
    /// `1e300` — 1e300 chromosomes, nothing the caller can meet, and far past the 9.0e15 the
    /// repeat-tract seed's own worst case reaches. So this test runs the two inputs that reach it:
    /// a legal-but-absurd total, and the ordinary case of an alternative entry sitting exactly on
    /// `MIN_ALT_CONCENTRATION`, which a fully invariant cohort produces.
    #[test]
    fn the_frequency_floor_keeps_every_entry_finite_at_the_edges() {
        // A fully invariant cohort's seed: the alternative concentration floored, the reference
        // ordinary. Reachable, and the row must be finite and ordered.
        let invariant = [1.0, crate::genetics::MIN_ALT_CONCENTRATION];
        let row = row_under(&PlugInWrightPrior, &invariant, 2, 0.0);
        assert!(row[0] > row[1] && row[1] > row[2], "{row:?}");

        // A total large enough that the floor is what keeps the row finite. `Concentration::new`
        // accepts this: both entries are finite and at or above `MIN_ALT_CONCENTRATION`.
        let absurd = [1e300, crate::genetics::MIN_ALT_CONCENTRATION];
        let row = row_under(&PlugInWrightPrior, &absurd, 2, 0.0);
        // `row_under` already refuses a non-finite entry, so reaching here is the assertion; what
        // this adds is that the floored entries are the floor's own value and not something the
        // arithmetic wandered to.
        let floored = 2.0 * PROBABILITY_FLOOR.ln();
        assert!(
            (row[2] - floored).abs() < 1.0,
            "{} against {floored}",
            row[2]
        );

        // And at four copies rather than two, where the floored frequency is multiplied by more.
        let row = row_under(&PlugInWrightPrior, &absurd, 4, 0.0);
        assert!(row.iter().all(|entry| entry.is_finite()), "{row:?}");
    }

    /// **The comparator allocates nothing**, like the default — it fills the caller's row and the
    /// caller's per-allele scratch and holds no buffer of its own. Pinned by the shape of the
    /// call rather than by a counter: `Cargo.toml` forbids `unsafe_code` crate-wide, so a
    /// counting allocator cannot be installed (recorded in `calling_loop.md`).
    ///
    /// What this checks instead is the observable consequence: the same buffers, reused across
    /// two calls at different inbreeding coefficients, give the two rows those coefficients imply
    /// rather than a mixture of them.
    #[test]
    fn the_caller_owns_every_buffer_and_they_survive_reuse() {
        let concentration = [1.0, TOMATO_DIVERSITY];
        let table = GenotypeTable::build(Ploidy::try_new(2).unwrap(), 2);
        let view = table.view();
        let mut scratch = vec![0.0; 2];
        let mut out = vec![LogProb(f64::NAN); view.genotype_count()];

        let mut first = Vec::new();
        for inbreeding in [0.0, 0.8, 0.0] {
            let mut row = PriorRow::new(
                Concentration::new(&concentration),
                view.genotype_allele_counts(),
                view.log_multinomial_coeffs(),
                view.homozygous_alleles(),
                &mut scratch,
                &mut out,
            );
            PlugInWrightPrior
                .fill_genotype_log_priors(&mut row, InbreedingF::try_new(inbreeding).unwrap());
            let taken: Vec<f64> = out.iter().map(|entry| entry.get()).collect();
            if first.is_empty() {
                first = taken;
            } else if inbreeding == 0.0 {
                // The third call must reproduce the first bit for bit: nothing carried over.
                assert_eq!(taken, first);
            } else {
                assert_ne!(taken, first);
            }
        }
    }
}
