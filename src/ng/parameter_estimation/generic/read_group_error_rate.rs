//! One error rate per read group, fitted from that group's own table.
//!
//! **The `ε` half of the coupled alternation** (`spec/parameter_prepass_generic.md` §5.1).
//! Each read group's rate is scanned over the error-rate ladder with the sample's genotype
//! frequencies held where the caller put them — one shared set across the groups, because
//! the frequencies are a property of the individual while the rates are a property of the
//! chemistry. Nothing here climbs the frequencies; the other half of the alternation does
//! that, on the whole-sample table, and hands the answer back for the next round
//! (Milestone E2).
//!
//! **Why a read group's table is scored as though that group were the sample.** The
//! read-group table enters a site once per group that covered it, at that group's own
//! depth and its own alternative count — so every cell in it belongs to one library and
//! the multi-library sum of `noise_model`'s expression collapses to a single term at share
//! one. [`SampleLibraryNoise::single`] is that statement, and it is where the group being
//! fitted and the rate being tried are tied together in one value rather than travelling
//! as a pair.
//!
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §5.1, spec §3.

use std::collections::BTreeMap;

use smallvec::SmallVec;

use crate::ng::parameter_estimation::fitting::ladder_scan::fit_by_fixed_frequency_scan;
use crate::ng::parameter_estimation::generic::histogram::{Cell, DepthAltHistogram};
use crate::ng::parameter_estimation::generic::noise_model::{
    SampleLibraryNoise, SubstitutionNoiseModel,
};
use crate::ng::types::{ErrorRate, LogProb, Ploidy, ReadGroupId};

/// What one read group's error-rate scan returned.
///
/// **Not an [`Estimate`](crate::ng::parameter_estimation::Estimate).** Everything here was
/// fitted from this group's own table, so there is no provenance to carry yet; the
/// fallback ladder that decides between a fitted rate, a borrowed one, a supplied one and
/// the default is Milestone E4's, and it is what attaches a
/// [`Provenance`](crate::ng::parameter_estimation::Provenance).
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct ReadGroupErrorRate {
    /// The winning rung's error rate.
    pub error_rate: ErrorRate,
    /// Where that rung sat on the ladder, counting from zero.
    ///
    /// **The coupled fit's stopping rule is stated in rungs** — it stops when every read
    /// group's winning rung is the one it had last iteration
    /// (`arch/parameter_prepass_generic.md` §5.2) — so this travels rather than being
    /// re-derived by searching the ladder for `error_rate`.
    pub rung: usize,
    /// The weighted log-likelihood of this group's cells at the winning rung and at the
    /// genotype frequencies the fit was handed. Comparable across the rungs of one call;
    /// across two calls only if the frequencies were the same.
    pub log_likelihood: LogProb,
    /// Whether the answer sat on the ladder's edge. **The one bit between a railed fit and
    /// a plausible-looking number**: a group whose true rate lies outside Phred 10–50 has
    /// its answer clamped to an endpoint and reported with the group's whole site count
    /// behind it.
    pub argmax_at_ladder_end: bool,
    /// How many sites this group's fit was made from, summed over every ploidy it covered.
    ///
    /// Carried because [`MIN_SITES_TO_FIT`](super::MIN_SITES_TO_FIT) is checked against it
    /// in Milestone E4, and the alternative is for that step to re-walk the tables to
    /// recover a number this one already summed.
    pub sites: u64,
}

/// Fit every read group's error rate from its own table, at the genotype frequencies
/// handed in.
///
/// `read_group_histograms` is the read-group table as the accumulator holds it — keyed by
/// `(read group, ploidy)`, one histogram each
/// ([`GenericAccumulators::read_group_histograms`](super::accumulators::GenericAccumulators::read_group_histograms)).
/// A group's cells are **every ploidy it covered, gathered into one scan**: one error rate
/// is fitted across all of them, because a haploid sex chromosome and the diploid
/// autosomes were prepared by the same chemistry, and each cell is scored against its own
/// genotype set.
///
/// `genotype_frequencies` carries one set per ploidy, each summing to one, in ascending
/// order of alternative copies. **One shared set across the groups** — that is what the
/// research harness's alternation does in both its arms, and it is what makes the two
/// blocks of the coupled fit two blocks rather than one per library.
///
/// `ladder` is the rungs to try, ordinarily [`error_rate_ladder`](super::error_rate_ladder).
/// Taken rather than built here so that a caller alternating twenty times builds it once,
/// and so that a test can hand in a short ladder whose right answer it can state.
///
/// **A read group holding no cells at all is left out of the answer** rather than fitted
/// or defaulted. It cannot arise through `add_locus`, which creates a histogram only when
/// a site enters it; if it does arise, there is nothing to fit, and Milestone E4's fallback
/// ladder is the one place that decides what a group with no fit gets — a borrowed rate or
/// the default — rather than this function inventing one.
///
/// # Panics
///
/// If `ladder` is empty, if `genotype_frequencies` has no entry for a ploidy some group
/// covered, if an entry is not as wide as that ploidy's genotype set, or if an entry is not
/// a probability vector. Also if a cell of a group's table lists a read group other than
/// that group's own, which would mean the read-group table and these parameters were built
/// from two different samples.
#[must_use]
pub fn fit_read_group_error_rates(
    read_group_histograms: &BTreeMap<(ReadGroupId, Ploidy), DepthAltHistogram<u64>>,
    genotype_frequencies: &BTreeMap<Ploidy, SmallVec<[f64; 3]>>,
    ladder: &[ErrorRate],
) -> BTreeMap<ReadGroupId, ReadGroupErrorRate> {
    assert!(!ladder.is_empty(), "a scan needs at least one rung to try");

    // The map is keyed `(group, ploidy)` and `BTreeMap` orders by the pair, so a group's
    // ploidies arrive together and in order — which is what lets one pass gather each
    // group's cells without a second index.
    let mut cells_of_group: BTreeMap<ReadGroupId, Vec<Cell>> = BTreeMap::new();
    for (&(read_group, ploidy), histogram) in read_group_histograms {
        cells_of_group
            .entry(read_group)
            .or_default()
            .extend(histogram.cells(ploidy));
    }

    let model = SubstitutionNoiseModel;
    let mut fitted = BTreeMap::new();
    for (read_group, cells) in cells_of_group {
        if cells.is_empty() {
            continue;
        }
        let sites = cells.iter().map(|cell| cell.sites).sum();

        // One rung, one parameter set, and the group travels **inside** it. Building the
        // ladder here rather than once outside the loop is what keeps the rate being tried
        // and the group it belongs to from being two collections indexed by position —
        // the fault `spec/parameter_prepass_generic.md` §1's scoring rule has no identity
        // against, since a rule with two libraries' rates swapped is still a probability.
        let noise_ladder: Vec<SampleLibraryNoise> = ladder
            .iter()
            .map(|&error_rate| SampleLibraryNoise::single(read_group, error_rate))
            .collect();

        let scanned =
            fit_by_fixed_frequency_scan(&model, &cells, &noise_ladder, genotype_frequencies);

        let libraries = scanned.noise.libraries();
        assert_eq!(
            libraries.len(),
            1,
            "read group {} was scanned against {} libraries, where its own table has one \
             library by construction",
            read_group.get(),
            libraries.len()
        );
        assert_eq!(
            libraries[0].read_group,
            read_group,
            "read group {} was fitted with a rate belonging to read group {}",
            read_group.get(),
            libraries[0].read_group.get()
        );

        fitted.insert(
            read_group,
            ReadGroupErrorRate {
                error_rate: libraries[0].error_rate,
                rung: scanned.rung,
                log_likelihood: scanned.log_likelihood,
                argmax_at_ladder_end: scanned.argmax_at_ladder_end,
                sites,
            },
        );
    }
    fitted
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
    use crate::ng::parameter_estimation::generic::error_rate_ladder;
    use crate::ng::parameter_estimation::generic::histogram::DepthAndAltReads;
    use crate::ng::types::Bp;

    /// Phred 30, which is rung 80 of the adopted ladder — the ladder starts at Phred 10
    /// and steps by 0.25, so rung `r` is Phred `10 + r/4`.
    const RUNG_AT_PHRED_30: usize = 80;
    /// Phred 26, rung 64 — a second library that is genuinely worse than the first, four
    /// Phred and sixteen rungs away, so a fit that answered one group from the other's
    /// table lands sixteen rungs out rather than one.
    const RUNG_AT_PHRED_26: usize = 64;

    fn ploidy(copies: u8) -> Ploidy {
        Ploidy::try_new(copies).expect("a positive copy number")
    }

    fn frequencies(entries: &[(u8, &[f64])]) -> BTreeMap<Ploidy, SmallVec<[f64; 3]>> {
        entries
            .iter()
            .map(|&(copies, set)| (ploidy(copies), SmallVec::from_slice(set)))
            .collect()
    }

    /// `p_j(ε)` — the same expression `noise_model` scores with, restated here rather
    /// than reached for, so that the fixture and the model are two statements of the rule
    /// and a fit that recovers its own rung is a claim about the **scan** rather than
    /// about the expression. Whether the expression itself is right is D2's four
    /// identities, not this file's.
    fn alternative_read_probability(alt_copies: u8, ploidy: Ploidy, error_rate: f64) -> f64 {
        let carried = f64::from(alt_copies) / f64::from(ploidy.get());
        carried * (1.0 - error_rate / 3.0) + (1.0 - carried) * error_rate
    }

    /// A table of `sites` sites, every one at `depth`, filled with **the counts an
    /// infinite genome would produce** at these parameters: cell `k` gets
    /// `sites · Σ_j π_j · Binomial(k; depth, p_j(ε))`, rounded.
    ///
    /// The same device the research harnesses use — replace the observed count with its
    /// probability under a known truth and the answer has no sampling noise in it, so a
    /// departure from the generating rung is bias rather than luck (research note §1).
    /// One depth, so the cell's mean depth is that depth exactly and the binning rule
    /// contributes nothing to the answer.
    fn table_generated_at(
        edges: &Arc<DepthBinEdges>,
        depth: u32,
        error_rate: f64,
        ploidy: Ploidy,
        genotype_frequencies: &[f64],
        sites: f64,
    ) -> DepthAltHistogram<u64> {
        let mut histogram = DepthAltHistogram::new(Arc::clone(edges));
        for alt_reads in 0..=depth {
            let mut probability = 0.0;
            for (alt_copies, &frequency) in genotype_frequencies.iter().enumerate() {
                let p = alternative_read_probability(alt_copies as u8, ploidy, error_rate);
                // The binomial term, built up from `k = 0` rather than through a
                // factorial, which at depth 20 keeps every intermediate inside `f64`.
                let mut term = (1.0 - p).powi(depth as i32);
                for step in 1..=alt_reads {
                    term *= f64::from(depth - step + 1) / f64::from(step) * p / (1.0 - p);
                }
                probability += frequency * term;
            }
            let in_cell = (sites * probability).round() as u64;
            for _ in 0..in_cell {
                histogram.add_site(DepthAndAltReads::new(depth, alt_reads), Bp(1));
            }
        }
        histogram
    }

    /// Tomato-like: heterozygous at 1.5 sites in a thousand, homozygous non-reference at
    /// 0.5, and everything else homozygous reference.
    const TRUTH: [f64; 3] = [0.998, 0.0015, 0.0005];
    const DEPTH: u32 = 20;
    const SITES: f64 = 200_000.0;

    /// **The claim E1 exists for, and the one that cannot be made from a single group:
    /// two libraries of one sample, four Phred apart, each recover their own rate.** A
    /// fit that answered either group from the other's table lands sixteen rungs out, and
    /// a fit that pooled the two tables lands between them — neither is a near miss.
    #[test]
    fn two_read_groups_at_different_rates_each_recover_their_own() {
        let edges = Arc::new(DepthBinEdges::new());
        let ladder = error_rate_ladder();
        let diploid = ploidy(2);

        let mut histograms = BTreeMap::new();
        for (group, rung) in [(1u32, RUNG_AT_PHRED_30), (2, RUNG_AT_PHRED_26)] {
            histograms.insert(
                (ReadGroupId(group), diploid),
                table_generated_at(&edges, DEPTH, ladder[rung].get(), diploid, &TRUTH, SITES),
            );
        }

        let fitted = fit_read_group_error_rates(&histograms, &frequencies(&[(2, &TRUTH)]), &ladder);

        assert_eq!(fitted.len(), 2);
        assert_eq!(fitted[&ReadGroupId(1)].rung, RUNG_AT_PHRED_30);
        assert_eq!(fitted[&ReadGroupId(2)].rung, RUNG_AT_PHRED_26);
        for group in [1, 2] {
            let fit = &fitted[&ReadGroupId(group)];
            assert!(
                !fit.argmax_at_ladder_end,
                "group {group} won at rung {}, which is interior",
                fit.rung
            );
            assert_eq!(fit.error_rate, ladder[fit.rung]);
        }
    }

    /// **The frequencies are held fixed, and holding them somewhere else moves the
    /// rate** — which is the whole reason E1 is the `ε` half of an alternation rather
    /// than a fit that stands alone. Told the sample is ten times as heterozygous as it
    /// is, the scan explains the same alternative reads as real variation and takes a
    /// **lower** error rate, which on this ladder is a **higher** rung.
    ///
    /// A scan that ignored the frequencies, or that climbed its own at every rung, would
    /// return the same rung twice.
    ///
    /// **Three reads a site, not twenty, and the depth is what makes this test able to
    /// fail at all.** The coupling between the two parameters is a low-coverage
    /// phenomenon: at twenty reads a heterozygote shows about ten alternative reads and a
    /// homozygous-reference site about one in fifty, so no plausible heterozygosity
    /// competes for the one-alternative-read cell and the same fixture returns rung 80
    /// under both frequency sets — measured, before this test was moved down to three
    /// reads. At three, a heterozygote shows one alternative read of three three times in
    /// eight, which is where the two explanations of that cell overlap and where the
    /// research note fits its coupled worlds.
    #[test]
    fn the_frequencies_handed_in_move_the_fitted_rate() {
        const SHALLOW: u32 = 3;
        let edges = Arc::new(DepthBinEdges::new());
        let ladder = error_rate_ladder();
        let diploid = ploidy(2);
        let histograms = BTreeMap::from([(
            (ReadGroupId(1), diploid),
            table_generated_at(
                &edges,
                SHALLOW,
                ladder[RUNG_AT_PHRED_30].get(),
                diploid,
                &TRUTH,
                SITES,
            ),
        )]);

        let at_the_truth =
            fit_read_group_error_rates(&histograms, &frequencies(&[(2, &TRUTH)]), &ladder);
        let told_it_is_variable = fit_read_group_error_rates(
            &histograms,
            &frequencies(&[(2, &[0.98, 0.015, 0.005])]),
            &ladder,
        );

        assert_eq!(
            at_the_truth[&ReadGroupId(1)].rung,
            RUNG_AT_PHRED_30,
            "at the truth's own frequencies the generating rung comes back"
        );
        assert!(
            told_it_is_variable[&ReadGroupId(1)].rung > RUNG_AT_PHRED_30,
            "at ten times the true heterozygosity the fit stayed at rung {}",
            told_it_is_variable[&ReadGroupId(1)].rung
        );
    }

    /// **One rate per read group across every ploidy it covered** — chemistry does not
    /// know about chromosomes. The haploid cells are generated at the same rate as the
    /// diploid ones, so both arms agree about the answer; what this pins is that the two
    /// are one scan, with the site count of both behind it, rather than two entries or one
    /// arm silently dropped.
    #[test]
    fn a_group_covering_two_ploidies_is_fitted_once_across_both() {
        let edges = Arc::new(DepthBinEdges::new());
        let ladder = error_rate_ladder();
        let rate = ladder[RUNG_AT_PHRED_30].get();
        let haploid_sites = 40_000.0;

        let histograms = BTreeMap::from([
            (
                (ReadGroupId(1), ploidy(2)),
                table_generated_at(&edges, DEPTH, rate, ploidy(2), &TRUTH, SITES),
            ),
            (
                (ReadGroupId(1), ploidy(1)),
                table_generated_at(
                    &edges,
                    DEPTH,
                    rate,
                    ploidy(1),
                    &[0.999, 0.001],
                    haploid_sites,
                ),
            ),
        ]);

        let fitted = fit_read_group_error_rates(
            &histograms,
            &frequencies(&[(1, &[0.999, 0.001]), (2, &TRUTH)]),
            &ladder,
        );

        assert_eq!(fitted.len(), 1, "one rate, not one per ploidy");
        let fit = &fitted[&ReadGroupId(1)];
        assert_eq!(fit.rung, RUNG_AT_PHRED_30);
        assert_eq!(
            fit.sites,
            histograms
                .values()
                .map(DepthAltHistogram::total_loci)
                .sum::<u64>(),
            "the site count spans both ploidies"
        );
    }

    /// **The rail flag, which is the only thing between a railed fit and a
    /// plausible-looking number.** A library ten times noisier than the ladder's worst
    /// rung has its answer clamped to that rung and reported with the group's whole site
    /// count behind it. The contrast that gives this teeth is the interior winner of
    /// `two_read_groups_at_different_rates_each_recover_their_own`.
    #[test]
    fn a_group_noisier_than_the_ladder_reaches_sets_the_rail_flag() {
        let edges = Arc::new(DepthBinEdges::new());
        let ladder = error_rate_ladder();
        let diploid = ploidy(2);
        // Phred 5 — the ladder's worst rung is Phred 10, so the truth is off its end.
        let histograms = BTreeMap::from([(
            (ReadGroupId(1), diploid),
            table_generated_at(&edges, DEPTH, 10f64.powf(-0.5), diploid, &TRUTH, SITES),
        )]);

        let fitted = fit_read_group_error_rates(&histograms, &frequencies(&[(2, &TRUTH)]), &ladder);

        let fit = &fitted[&ReadGroupId(1)];
        assert_eq!(fit.rung, 0, "clamped to the noisiest rung the ladder has");
        assert!(fit.argmax_at_ladder_end);
    }

    /// A read group whose table holds no site at all is left out rather than fitted or
    /// defaulted: there is nothing to fit, and E4's fallback ladder is the one place that
    /// decides what such a group gets.
    #[test]
    fn a_read_group_with_no_cells_is_left_out() {
        let edges = Arc::new(DepthBinEdges::new());
        let ladder = error_rate_ladder();
        let diploid = ploidy(2);
        let histograms = BTreeMap::from([
            (
                (ReadGroupId(1), diploid),
                table_generated_at(
                    &edges,
                    DEPTH,
                    ladder[RUNG_AT_PHRED_30].get(),
                    diploid,
                    &TRUTH,
                    SITES,
                ),
            ),
            (
                (ReadGroupId(2), diploid),
                DepthAltHistogram::new(Arc::clone(&edges)),
            ),
        ]);

        let fitted = fit_read_group_error_rates(&histograms, &frequencies(&[(2, &TRUTH)]), &ladder);

        assert_eq!(fitted.keys().copied().collect::<Vec<_>>(), [ReadGroupId(1)]);
    }

    #[test]
    #[should_panic(expected = "a scan needs at least one rung to try")]
    fn an_empty_ladder_is_refused() {
        let edges = Arc::new(DepthBinEdges::new());
        let diploid = ploidy(2);
        let histograms = BTreeMap::from([(
            (ReadGroupId(1), diploid),
            table_generated_at(&edges, DEPTH, 0.001, diploid, &TRUTH, 1_000.0),
        )]);

        let _ = fit_read_group_error_rates(&histograms, &frequencies(&[(2, &TRUTH)]), &[]);
    }
}
