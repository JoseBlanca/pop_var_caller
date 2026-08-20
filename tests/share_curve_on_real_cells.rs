//! The two slippage shares' curves, fitted against the strata two real cohorts produced.
//!
//! **The unit tests in `share_curve.rs` draw strata on a known shape and check the fit finds it
//! back. This one has no known answer to find** — it runs the shape choice over the per-stratum
//! tables the walk wrote for the tomato cohort and for HG002, and asserts what
//! `doc/devel/ng/reports/str_slippage_share_families_2026-08-20.md` reports: that no single shape
//! serves every motif period, and that all three earn their place.
//!
//! **Both tables were produced with every stratum fitted from its own tracts and nothing copied
//! from a neighbour.** Reading a copied share as a measurement would fit a curve to another
//! curve's output, which is the circularity the design forbids. On these two tables nothing was
//! copied: no tomato stratum ever cleared the copy rule's floor, and the HG002 table is the one
//! taken before that rule existed.
//!
//! Fixtures: `tests/data/slippage_cells/`, at the census's corrected ±8 recording window.

use std::collections::BTreeMap;
use std::path::Path;

use pop_var_caller::ng::parameter_estimation::joint::share_curve::{
    FittedShare, ShareCurveConfig, ShareCurveSource, ShareShape, choose_share_shape,
};

/// Which of the two shares a column holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Share {
    /// Of the reads that slipped, the share that came back shorter.
    DirectionSplit,
    /// How much rarer a two-unit slip is than a one-unit slip.
    FallOff,
}

impl Share {
    fn column(self) -> &'static str {
        match self {
            Self::DirectionSplit => "shorter_share",
            Self::FallOff => "fall_off",
        }
    }
}

/// One motif period's contributing strata for one share, read out of a cell table.
///
/// **Only rows the walk marked as fitted are here**, and the slipped-read count is the stratum's
/// own fitted level times the reads that crossed it — the count the two shares rest on.
fn strata_by_period(path: &Path, share: Share) -> BTreeMap<u8, Vec<FittedShare>> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|error| panic!("{} is readable: {error}", path.display()));
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().expect("a header row").split(',').collect();
    let column = |name: &str| {
        header
            .iter()
            .position(|field| *field == name)
            .unwrap_or_else(|| panic!("the cell table has a {name} column"))
    };
    let (period, repeats, reads, fitted, level, wanted) = (
        column("period"),
        column("repeats"),
        column("spanning_reads"),
        column("fitted"),
        column("level"),
        column(share.column()),
    );

    let mut by_period: BTreeMap<u8, Vec<FittedShare>> = BTreeMap::new();
    for line in lines.filter(|line| !line.trim().is_empty()) {
        let field: Vec<&str> = line.split(',').collect();
        if field[fitted] != "1" {
            continue;
        }
        let level: f64 = field[level].parse().expect("a fitted level");
        let spanning: f64 = field[reads].parse().expect("a spanning-read count");
        by_period
            .entry(field[period].parse().expect("a motif period"))
            .or_default()
            .push(FittedShare {
                repeats: field[repeats].parse().expect("a repeat count"),
                share: field[wanted].parse().expect("a fitted share"),
                slipped_reads: level * spanning,
            });
    }
    for strata in by_period.values_mut() {
        strata.sort_by_key(|stratum| stratum.repeats);
    }
    by_period
}

fn fixture(name: &str, share: Share) -> BTreeMap<u8, Vec<FittedShare>> {
    strata_by_period(
        &Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/slippage_cells")
            .join(name),
        share,
    )
}

/// Every period of both cohorts that has enough fitted strata to choose a shape between.
///
/// Ten in all: four motif periods on HG002, one on tomato, two shares each.
fn every_case() -> Vec<(&'static str, u8, Share, usize, ShareShape, f64)> {
    let mut cases = Vec::new();
    for name in ["tomato_63_accessions_8mb.csv", "hg002_300x_tier.csv"] {
        for share in [Share::DirectionSplit, Share::FallOff] {
            for (period, strata) in fixture(name, share) {
                let Some(curve) = choose_share_shape(&strata, &ShareCurveConfig::default()) else {
                    continue;
                };
                cases.push((
                    name,
                    period,
                    share,
                    strata.len(),
                    curve.shape,
                    curve.held_out_error,
                ));
            }
        }
    }
    cases
}

/// **The finding the three shapes exist for: no one of them serves every period.** Drop any of
/// them and some period of some cohort is fitted with a shape its own strata reject.
#[test]
fn all_three_shapes_win_somewhere_across_the_two_cohorts() {
    let cases = every_case();
    assert_eq!(cases.len(), 10, "ten cases can be scored: {cases:#?}");
    for shape in ShareShape::ALL {
        assert!(
            cases.iter().any(|(_, _, _, _, won, _)| *won == shape),
            "{shape:?} never wins, so it could be dropped: {cases:#?}"
        );
    }
}

/// Which shape each period's strata choose, and how far that shape lands from a stratum it never
/// saw.
///
/// **These are asserted rather than left to the prose** so the report's table cannot go stale
/// silently. The errors are in logit units — the units the blend weighs a curve against a
/// stratum's own answer in — and they are pinned from this fit rather than predicted from
/// anywhere else.
///
/// **Two of the ten are near ties and the tie rule decides them**, which is worth knowing before
/// reading a shape as a fact about the chemistry: HG002's homopolymer direction split and its
/// dinucleotide fall-off are each within a few hundredths of a logit unit of the next shape up,
/// and the simpler shape keeps them.
#[test]
fn each_period_chooses_the_shape_the_report_records() {
    use Share::{DirectionSplit as Split, FallOff};

    let expected = [
        (
            "tomato_63_accessions_8mb.csv",
            1_u8,
            Split,
            5_usize,
            ShareShape::Sloping,
            0.240,
        ),
        (
            "tomato_63_accessions_8mb.csv",
            1,
            FallOff,
            5,
            ShareShape::Flat,
            0.321,
        ),
        ("hg002_300x_tier.csv", 1, Split, 23, ShareShape::Flat, 0.165),
        (
            "hg002_300x_tier.csv",
            1,
            FallOff,
            23,
            ShareShape::Sloping,
            0.239,
        ),
        (
            "hg002_300x_tier.csv",
            2,
            Split,
            20,
            ShareShape::Turning,
            0.516,
        ),
        (
            "hg002_300x_tier.csv",
            2,
            FallOff,
            20,
            ShareShape::Flat,
            0.527,
        ),
        (
            "hg002_300x_tier.csv",
            3,
            Split,
            4,
            ShareShape::Turning,
            0.210,
        ),
        (
            "hg002_300x_tier.csv",
            3,
            FallOff,
            4,
            ShareShape::Flat,
            0.813,
        ),
        ("hg002_300x_tier.csv", 4, Split, 7, ShareShape::Flat, 0.789),
        (
            "hg002_300x_tier.csv",
            4,
            FallOff,
            7,
            ShareShape::Flat,
            0.840,
        ),
    ];

    let cases = every_case();
    for (name, period, share, strata, shape, error) in expected {
        let found = cases
            .iter()
            .find(|(at, at_period, at_share, ..)| {
                *at == name && *at_period == period && *at_share == share
            })
            .unwrap_or_else(|| panic!("{name} period {period} {share:?} should be scorable"));
        assert_eq!(found.3, strata, "{name} period {period} {share:?} strata");
        assert_eq!(found.4, shape, "{name} period {period} {share:?} shape");
        assert!(
            (found.5 - error).abs() < 0.005,
            "{name} period {period} {share:?} predicts a held-out stratum to {:.3} logit units, \
             expected {error:.3}",
            found.5
        );
    }
}

/// **What a curve is worth against the strata it is meant to serve.** A well-measured stratum
/// holds its own share far more precisely than any period's curve predicts one — so the curve is
/// for the strata that have no answer of their own, and a stratum that does have one keeps most
/// of its weight.
#[test]
fn a_well_measured_stratum_holds_its_share_more_precisely_than_its_curve() {
    let strata = fixture("hg002_300x_tier.csv", Share::DirectionSplit);
    let homopolymers = &strata[&1];
    let curve = choose_share_shape(homopolymers, &ShareCurveConfig::default())
        .expect("twenty-three strata");

    let mut own: Vec<f64> = homopolymers
        .iter()
        .map(|stratum| stratum.logit_standard_error())
        .collect();
    own.sort_by(f64::total_cmp);
    let median_own = own[own.len() / 2];

    assert!(
        median_own * 4.0 < curve.held_out_error,
        "the median stratum holds its split to {median_own:.3} logit units and the curve \
         predicts one to {:.3}",
        curve.held_out_error
    );
}

/// Every period that can be scored comes back with a curve fitted from its own strata, and every
/// share it reports is a proportion at every repeat count a tract can have.
#[test]
fn every_scored_period_gives_shares_that_are_proportions() {
    for (name, period, share, ..) in every_case() {
        let strata = fixture(name, share);
        let curve = choose_share_shape(&strata[&period], &ShareCurveConfig::default())
            .expect("this period was scorable");
        assert_eq!(curve.source, ShareCurveSource::ThisPeriod);
        for repeats in 1..=200 {
            let value = curve.share_at(repeats);
            assert!(
                (0.0..=1.0).contains(&value) && value.is_finite(),
                "{name} period {period} {share:?} gave {value} at {repeats} repeats"
            );
        }
    }
}

/// **The one curve drawn through only four strata**, and the reason it is allowed. HG002's
/// trinucleotide direction split swings 2.7-fold across four repeat counts, and a shape that
/// turns is the only one of the three that follows it — with four strata and three coefficients
/// every leave-one-out fit passes exactly through the three that remain, so the score cannot
/// catch a bad shape here, and what stands in for that is how well it predicts against the other
/// nine cases.
#[test]
fn the_four_stratum_turning_curve_predicts_better_than_most_richer_ones() {
    let cases = every_case();
    let (.., strata, shape, error) = cases
        .iter()
        .find(|(name, period, share, ..)| {
            *name == "hg002_300x_tier.csv" && *period == 3 && *share == Share::DirectionSplit
        })
        .copied()
        .expect("HG002's trinucleotide direction split is scorable");
    assert_eq!(strata, 4);
    assert_eq!(shape, ShareShape::Turning);

    let worse = cases.iter().filter(|(.., other)| *other > error).count();
    assert_eq!(
        worse, 8,
        "eight of the other nine cases predict a held-out stratum worse than {error:.3}"
    );
}
