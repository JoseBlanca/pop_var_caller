use super::*;
use crate::neighbour_evidence::refine::MISMATCH;
use crate::neighbour_evidence::templates::{MAX_TEMPLATES, choose_templates};
use crate::neighbour_evidence::test_cohort::{CohortShape, TestCohort, test_cohort};

/// **`F` from the parameters file is clamped** to [0.01, 0.99], and a missing one (NaN) to the
/// lower bound.
#[test]
fn the_starting_share_is_clamped() {
    assert_eq!(
        SampleFit::new(1.0).parameters(MISMATCH).autozygous_share,
        0.99
    );
    assert_eq!(
        SampleFit::new(0.0).parameters(MISMATCH).autozygous_share,
        0.01
    );
    assert_eq!(
        SampleFit::new(f64::NAN)
            .parameters(MISMATCH)
            .autozygous_share,
        0.01
    );
    assert_eq!(
        SampleFit::new(0.7).parameters(MISMATCH).autozygous_share,
        0.7
    );
}

/// The walk sites of sample `sample` over the whole cohort, its templates and the background, as
/// the window solve builds them; returned as owned buffers and the template count.
/// One walk site's position, likelihoods and template shares, owned.
type OwnedSite = (u64, Vec<f64>, Vec<f64>);

fn sample_sites(cohort: &TestCohort, sample: usize) -> (Vec<OwnedSite>, usize) {
    let templates =
        choose_templates(&cohort.records, sample, MAX_TEMPLATES).expect("informative records");
    let sites = cohort
        .records
        .iter()
        .map(|record| {
            let mut shares: Vec<f64> = Vec::new();
            for &template in &templates {
                shares.extend(record.template_shares(template));
            }
            shares.extend(record.background_shares());
            (
                record.position,
                record.likelihoods(sample).collect(),
                shares,
            )
        })
        .collect();
    (sites, templates.len() + 1)
}

fn walk_sites(owned: &[OwnedSite]) -> Vec<WalkSite<'_>> {
    owned
        .iter()
        .map(|(position, likelihoods, shares)| WalkSite {
            position: *position,
            allele_count: 2,
            likelihoods,
            template_shares: shares,
        })
        .collect()
}

/// **The fit freezes once its windows hold [`FIT_RECORDS`] informative records**, and a frozen
/// fit no longer moves.
#[test]
fn the_fit_freezes_after_enough_informative_records() {
    let cohort = test_cohort(CohortShape {
        samples: 12,
        founders: 3,
        records: 200,
        spacing: 1_000,
        segment_bases: 100_000.0,
        inbred: true,
        depth: 3.0,
        read_error: 0.01,
        seed: 3,
        silent_from: None,
    });
    let (owned, templates) = sample_sites(&cohort, 0);
    let sites = walk_sites(&owned);
    let mut fit = SampleFit::new(0.9);
    fit.fit_window(&sites, templates, FIT_RECORDS / 2, MISMATCH);
    assert!(!fit.is_frozen());
    let frozen = fit.fit_window(&sites, templates, FIT_RECORDS / 2, MISMATCH);
    assert!(fit.is_frozen());
    let after = fit.fit_window(&sites[..50], templates, 50, MISMATCH);
    assert_eq!(frozen, after);
}

/// **An inbred sample's `F` is re-estimated high, and an outbred one's low**, from the same
/// starting value: on a drawn cohort of 20 samples at 3 reads a record, wholly inbred against
/// wholly outbred.
#[test]
fn the_autozygous_share_follows_the_sample() {
    let shape = CohortShape {
        samples: 20,
        founders: 4,
        records: 2_000,
        spacing: 500,
        segment_bases: 300_000.0,
        inbred: true,
        depth: 3.0,
        read_error: 0.01,
        seed: 11,
        silent_from: None,
    };
    let share = |inbred: bool| {
        let cohort = test_cohort(CohortShape { inbred, ..shape });
        let (owned, templates) = sample_sites(&cohort, 0);
        let mut fit = SampleFit::new(0.5);
        fit.fit_window(&walk_sites(&owned), templates, owned.len(), MISMATCH)
            .autozygous_share
    };
    let inbred = share(true);
    let outbred = share(false);
    assert!(inbred > 0.9, "inbred sample: F {inbred}");
    assert!(outbred < 0.2, "outbred sample: F {outbred}");
}
