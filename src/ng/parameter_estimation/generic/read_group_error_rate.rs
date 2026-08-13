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
//! Design: `doc/devel/ng/spec/parameter_prepass_generic.md` §3 for what this table holds
//! and why only `ε` comes out of it, and §5.1 for the alternation. **Not §3's second
//! paragraph**, which describes the procedure this half replaced — "step through the error
//! rate, and at each step climb to the genotype frequencies that fit best" is
//! `fit_by_profile_scan`, and the owner's decision of 2026-08-07 is that the `ε` step holds
//! the frequencies fixed instead. `arch/parameter_prepass_generic.md` §5.1 and §5.2 carry
//! the same stale reading; that is recorded for the owner rather than edited here.

use std::collections::BTreeMap;

use crate::ng::parameter_estimation::fitting::ladder_scan::fit_by_fixed_frequency_scan;
use crate::ng::parameter_estimation::generic::SiteNoise;
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
pub struct ReadGroupErrorRateFit {
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
    ///
    /// **Sites, and `arch/parameter_prepass_generic.md` §2.4 says an error rate's
    /// `Estimate::observations` is *reads*.** The two differ by the mean depth, so
    /// whichever E4 puts on the estimate has to be the one the field's doc claims.
    /// Recorded here rather than converted, because a per-read count is a division this
    /// step has no reason to do and E4 has the histograms it would need.
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
/// `site_noise` is the sample's second class of site, when it has one, and **every rung is
/// scored with it.** Leaving it out was a defect rather than a simplification: a candidate
/// clean rate scored under the one-class rule is being asked to explain a table whose tail
/// belongs to the other class, so the scan returns the tail-inflated rate no matter what pair
/// sits beside it — and since this scan is where the coupled fit's rate comes from, the rate
/// then never moves. Measured on a world generated at HG002's own parameters it came back
/// **three rungs high**, the same rung the one-class fit chose, on that world and on all five
/// real alignments, with the whole fit scoring 351 nats below the truth
/// (`reports/implementations/ng_noise_model_extension_n5_2026-08-10.md`).
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
/// a probability vector.
///
/// And, two frames down in [`SampleLibraryNoise`], if a cell lists a read group other than
/// the one being fitted. **That reaches only the attributed arm**, and the read-group table
/// is built entirely of pooled keys, so on any table the accumulator produces it cannot
/// fire — which is exactly why the read group inside the noise parameters is inert here and
/// why the ladder is built inside the loop below rather than once outside it.
#[must_use]
pub fn fit_read_group_error_rates(
    read_group_histograms: &BTreeMap<(ReadGroupId, Ploidy), DepthAltHistogram<u64>>,
    genotype_frequencies: &BTreeMap<Ploidy, Vec<f64>>,
    ladder: &[ErrorRate],
    site_noise: Option<SiteNoise>,
) -> BTreeMap<ReadGroupId, ReadGroupErrorRateFit> {
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
        //
        // **And that construction is the whole of the guarantee, because the label is
        // inert on this path.** Every cell of the read-group table carries a pooled key,
        // and the pooled branch of the scoring rule reads a library's share and its rate
        // and never its read group — so a ladder hoisted out of this loop and built from
        // the first group's id fits every group to the right answer, and the two checks
        // below are what catch it. They compare against the same `read_group` the ladder
        // was built from, which is the point: they hold by construction while the
        // construction stays here, and stop holding the moment it moves.
        let noise_ladder: Vec<SampleLibraryNoise> = ladder
            .iter()
            .map(|&error_rate| {
                SampleLibraryNoise::single_with_site_noise(read_group, error_rate, site_noise)
            })
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
            ReadGroupErrorRateFit {
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
    use crate::ng::parameter_estimation::generic::expected_counts::table_generated_at;
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

    fn frequencies(entries: &[(u8, &[f64])]) -> BTreeMap<Ploidy, Vec<f64>> {
        entries
            .iter()
            .map(|&(copies, set)| (ploidy(copies), set.to_vec()))
            .collect()
    }

    /// Tomato-like: heterozygous at 1.5 sites in a thousand, homozygous non-reference at
    /// 0.5, and everything else homozygous reference.
    const TRUTH: [f64; 3] = [0.998, 0.0015, 0.0005];
    const DEPTH: u32 = 20;
    const SITES: f64 = 200_000.0;

    /// **The claim E1 exists for, and the one that cannot be made from a single group:
    /// two libraries of one sample, four Phred apart, each recover their own rate.** A
    /// fit that answered either group from the other's table lands sixteen rungs out, and
    /// a fit that pooled the two tables lands between them — neither is a near miss. The
    /// pooled arm is asserted rather than argued, because "between them" is a claim about
    /// this fixture and not a general fact.
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

        let fitted =
            fit_read_group_error_rates(&histograms, &frequencies(&[(2, &TRUTH)]), &ladder, None);

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

        // What a fit that pooled the two tables would return: one group's worth of sites
        // from each rate, entered under a single read group. Rung 70 is between 64 and 80
        // and is no library's rate — the failure has no symptom other than this number.
        let mut pooled = table_generated_at(
            &edges,
            DEPTH,
            ladder[RUNG_AT_PHRED_30].get(),
            diploid,
            &TRUTH,
            SITES,
        );
        pooled.merge(&table_generated_at(
            &edges,
            DEPTH,
            ladder[RUNG_AT_PHRED_26].get(),
            diploid,
            &TRUTH,
            SITES,
        ));
        let pooled = fit_read_group_error_rates(
            &BTreeMap::from([((ReadGroupId(1), diploid), pooled)]),
            &frequencies(&[(2, &TRUTH)]),
            &ladder,
            None,
        );
        assert_eq!(
            pooled[&ReadGroupId(1)].rung,
            70,
            "pooling the two libraries lands between their rungs, at neither one's rate"
        );
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
    /// **Six reads a site, not twenty, and the depth is what makes this test able to fail
    /// at all.** The coupling between the two parameters is a low-coverage phenomenon: at
    /// twenty reads a heterozygote shows about ten alternative reads (10.007) and a
    /// homozygous-reference site about one in fifty (one in 50.5), so no plausible
    /// heterozygosity competes for the one-alternative-read cell. Measured on this exact
    /// fixture, both frequency sets return rung 80 at twenty reads **and at ten**; the
    /// answer only moves below that.
    ///
    /// **Six and not three, and that is the second thing this fixture had to get
    /// right.** At three reads the second arm lands on rung **160**, the top of the
    /// ladder — a railed answer, which satisfies "higher than rung 80" without being an
    /// argmax of anything, so the assertion would also pass a defect that railed high for
    /// any frequency set but the fixture's own. At six the move is 80 → 84 and interior,
    /// which is why both halves are asserted below.
    #[test]
    fn the_frequencies_handed_in_move_the_fitted_rate() {
        const SHALLOW: u32 = 6;
        /// Where the second arm lands: four rungs, one Phred, below the true rate.
        const RUNG_TOLD_TEN_TIMES_AS_VARIABLE: usize = 84;
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
            fit_read_group_error_rates(&histograms, &frequencies(&[(2, &TRUTH)]), &ladder, None);
        let told_it_is_variable = fit_read_group_error_rates(
            &histograms,
            &frequencies(&[(2, &[0.98, 0.015, 0.005])]),
            &ladder,
            None,
        );

        assert_eq!(
            at_the_truth[&ReadGroupId(1)].rung,
            RUNG_AT_PHRED_30,
            "at the truth's own frequencies the generating rung comes back"
        );
        let moved = &told_it_is_variable[&ReadGroupId(1)];
        assert_eq!(
            moved.rung, RUNG_TOLD_TEN_TIMES_AS_VARIABLE,
            "at ten times the true heterozygosity the fit landed on rung {}",
            moved.rung
        );
        assert!(
            !moved.argmax_at_ladder_end,
            "the moved answer must be an argmax and not the ladder's end"
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
            None,
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
    /// plausible-looking number.** A library at Phred 5 — an error rate of 0.316, **3.2
    /// times** the 0.1 of the ladder's worst rung and five Phred past its end — has its
    /// answer clamped to that rung and reported with the group's whole site count behind
    /// it. The contrast that gives this teeth is the interior winner of
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

        let fitted =
            fit_read_group_error_rates(&histograms, &frequencies(&[(2, &TRUTH)]), &ladder, None);

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

        let fitted =
            fit_read_group_error_rates(&histograms, &frequencies(&[(2, &TRUTH)]), &ladder, None);

        assert_eq!(fitted.keys().copied().collect::<Vec<_>>(), [ReadGroupId(1)]);
    }

    /// **A group too thin to be worth fitting is still fitted here.**
    /// [`MIN_SITES_TO_FIT`](super::MIN_SITES_TO_FIT) is Milestone E4's gate, and E4 needs
    /// both the fit and the site count to decide against them — a group dropped at this
    /// step is one E4 cannot tell from a group that never had a table, and the two get
    /// different rungs of the fallback ladder.
    ///
    /// **Five hundred sites, against a floor of ten thousand**, because every other
    /// fixture in this file holds 40,000 or 200,000: a gate run a milestone early would be
    /// invisible to all of them.
    #[test]
    fn a_group_far_below_the_fitting_floor_is_still_returned_with_its_site_count() {
        let edges = Arc::new(DepthBinEdges::new());
        let ladder = error_rate_ladder();
        let diploid = ploidy(2);
        let thin = 500.0;
        assert!(
            thin < super::super::MIN_SITES_TO_FIT as f64,
            "the fixture has to sit below the floor for this test to say anything"
        );
        let histograms = BTreeMap::from([(
            (ReadGroupId(1), diploid),
            table_generated_at(
                &edges,
                DEPTH,
                ladder[RUNG_AT_PHRED_30].get(),
                diploid,
                &TRUTH,
                thin,
            ),
        )]);

        let fitted =
            fit_read_group_error_rates(&histograms, &frequencies(&[(2, &TRUTH)]), &ladder, None);

        let fit = fitted
            .get(&ReadGroupId(1))
            .expect("a thin group is fitted here and gated in E4");
        assert_eq!(
            fit.sites,
            histograms[&(ReadGroupId(1), diploid)].total_loci()
        );
    }

    /// **A cell is scored against the genotypes of the ploidy its table is keyed by.** A
    /// haploid-only group is handed haploid frequencies and nothing else, so a fit that
    /// labelled its cells with any other ploidy has no frequencies to look up and panics
    /// rather than answering.
    ///
    /// The two-ploidy test above cannot say this: it hands over both sets, so a fit that
    /// labelled every cell diploid finds a set waiting for it.
    #[test]
    fn a_haploid_group_is_scored_against_haploid_genotypes() {
        let edges = Arc::new(DepthBinEdges::new());
        let ladder = error_rate_ladder();
        let haploid = ploidy(1);
        let histograms = BTreeMap::from([(
            (ReadGroupId(1), haploid),
            table_generated_at(
                &edges,
                DEPTH,
                ladder[RUNG_AT_PHRED_30].get(),
                haploid,
                &[0.999, 0.001],
                SITES,
            ),
        )]);

        let fitted = fit_read_group_error_rates(
            &histograms,
            &frequencies(&[(1, &[0.999, 0.001])]),
            &ladder,
            None,
        );

        assert_eq!(fitted[&ReadGroupId(1)].rung, RUNG_AT_PHRED_30);
    }

    /// **The rate a group is scored against must be labelled with that group**, and this
    /// is the only test that can say so — with a caveat worth stating rather than hiding.
    ///
    /// The read-group table's cells are all **pooled**, and the pooled branch of the
    /// scoring rule reads a library's share and its rate and never its label. So on every
    /// table `add_locus` can build, the read group inside
    /// [`SampleLibraryNoise`] is **inert**: a ladder hoisted out of the per-group loop and
    /// built from the first group's id fits every group correctly. The two `assert_eq!`s
    /// in the fit are what catch that, and they can only catch it because the construction
    /// sits inside the loop.
    ///
    /// The **attributed** branch does read the label, and refuses a cell naming a library
    /// it does not have. So this test puts one attributed site in each group's table to
    /// give the pairing an observable consequence. `add_locus` never produces such a cell
    /// in the read-group table — that arm belongs to the windowed one — so what is pinned
    /// here is the function's contract rather than a reachable input.
    #[test]
    fn each_group_is_scored_against_noise_that_names_it() {
        let edges = Arc::new(DepthBinEdges::new());
        let ladder = error_rate_ladder();
        let diploid = ploidy(2);

        let mut histograms = BTreeMap::new();
        for (group, rung) in [(1u32, RUNG_AT_PHRED_30), (2, RUNG_AT_PHRED_26)] {
            let mut table =
                table_generated_at(&edges, DEPTH, ladder[rung].get(), diploid, &TRUTH, SITES);
            table.add_attributed_site(
                DepthAndAltReads::new(DEPTH, 1),
                &[(ReadGroupId(group), 1)],
                Bp(1),
            );
            histograms.insert((ReadGroupId(group), diploid), table);
        }

        let fitted =
            fit_read_group_error_rates(&histograms, &frequencies(&[(2, &TRUTH)]), &ladder, None);

        assert_eq!(fitted[&ReadGroupId(1)].rung, RUNG_AT_PHRED_30);
        assert_eq!(fitted[&ReadGroupId(2)].rung, RUNG_AT_PHRED_26);
    }

    /// **The score reported is the winning rung's own score**, which is the field the
    /// coupled fit of E2 picks its best-scoring iterate by — so a constant there would be
    /// read as "no iterate improved on the first".
    ///
    /// Two one-rung ladders make it checkable from outside: the rung the table was
    /// generated at must score above a rung sixteen away, and the full scan must report
    /// exactly the number its winning rung scored on its own.
    #[test]
    fn the_reported_score_is_the_winning_rungs_score() {
        let edges = Arc::new(DepthBinEdges::new());
        let ladder = error_rate_ladder();
        let diploid = ploidy(2);
        let histograms = BTreeMap::from([(
            (ReadGroupId(1), diploid),
            table_generated_at(
                &edges,
                DEPTH,
                ladder[RUNG_AT_PHRED_30].get(),
                diploid,
                &TRUTH,
                SITES,
            ),
        )]);
        let frequencies = frequencies(&[(2, &TRUTH)]);
        let score_at = |rung: usize| {
            fit_read_group_error_rates(&histograms, &frequencies, &ladder[rung..=rung], None)
                [&ReadGroupId(1)]
                .log_likelihood
                .get()
        };

        let at_the_truth = score_at(RUNG_AT_PHRED_30);
        assert!(
            at_the_truth < 0.0 && at_the_truth.is_finite(),
            "a weighted log-likelihood over {SITES} sites is a negative finite number, \
             not {at_the_truth}"
        );
        assert!(
            at_the_truth > score_at(RUNG_AT_PHRED_26),
            "the generating rung scored {at_the_truth} and a rung sixteen away scored {}",
            score_at(RUNG_AT_PHRED_26)
        );

        let whole_ladder = fit_read_group_error_rates(&histograms, &frequencies, &ladder, None);
        assert_eq!(
            whole_ladder[&ReadGroupId(1)].log_likelihood.get(),
            at_the_truth,
            "the scan reported a score its own winning rung does not have"
        );
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

        let _ = fit_read_group_error_rates(&histograms, &frequencies(&[(2, &TRUTH)]), &[], None);
    }
}
