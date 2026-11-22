use super::*;

fn biallelic(phred: &[u32; 3], posteriors: &[f64; 3]) -> Option<SampleColumns<'static>> {
    // Leaked so the columns can be built inline; a test's few bytes.
    Some(SampleColumns {
        phred_likelihoods: Box::leak(Box::new(*phred)),
        genotype_posteriors: Box::leak(Box::new(*posteriors)),
    })
}

/// **A record reads each sample's column as §3.2 says**: likelihoods `10^(−PL/10)`, allele shares
/// `Σ_g GP_g · copies(g) / 2`, a flat `PL` as a sample whose reads said nothing, a `.` column as
/// absent — likelihood one and equal shares.
#[test]
fn a_record_reads_likelihoods_and_shares_off_pl_and_gp() {
    let record = SiteRecord::new(
        100,
        &[0.25],
        &[
            biallelic(&[0, 10, 20], &[0.0, 1.0, 0.0]),
            biallelic(&[0, 0, 0], &[0.6, 0.3, 0.1]),
            None,
        ],
    );
    let close = |got: Vec<f64>, expected: &[f64]| {
        got.iter()
            .zip(expected)
            .all(|(got, expected)| (got - expected).abs() < 1e-6)
    };
    assert_eq!(record.row(0), SampleRow::Informative);
    assert_eq!(record.row(1), SampleRow::Flat);
    assert_eq!(record.row(2), SampleRow::Absent);
    assert!(close(record.likelihoods(0).collect(), &[1.0, 0.1, 0.01]));
    assert!(close(record.template_shares(0).collect(), &[0.5, 0.5]));
    assert!(close(record.template_shares(1).collect(), &[0.75, 0.25]));
    assert!(close(record.likelihoods(2).collect(), &[1.0, 1.0, 1.0]));
    assert!(close(record.template_shares(2).collect(), &[0.5, 0.5]));
    assert!(close(record.background_shares().collect(), &[0.75, 0.25]));
    assert!(record.is_site());
}

/// **A record is a site only with two samples' rows**: one to re-call and one to copy (§5).
#[test]
fn a_record_with_one_row_is_no_site() {
    let record = SiteRecord::new(
        100,
        &[0.25],
        &[biallelic(&[0, 10, 20], &[0.9, 0.1, 0.0]), None, None],
    );
    assert!(!record.is_site());
}

/// **Rounded frequencies that sum past one** leave the reference nothing and are rescaled.
#[test]
fn frequencies_past_one_are_rescaled() {
    let mut shares = Vec::new();
    background_shares(&[0.6, 0.6], &mut shares);
    assert_eq!(shares, vec![0.0, 0.5, 0.5]);
}

/// **A likelihood held as its `PL` above the lowest reads back as the `f32` the expression gives**,
/// bit for bit, at every `u16`, and a `PL` past `u16::MAX` reads as the zero it is.
#[test]
fn a_likelihood_held_in_two_bytes_is_the_f32_it_stands_for() {
    for phred in 0..=u16::MAX {
        let expected = crate::float::powf(10.0, -f64::from(phred) / 10.0) as f32;
        assert_eq!(
            likelihood_of(phred).to_bits(),
            expected.to_bits(),
            "PL {phred}"
        );
    }
    let far = crate::float::powf(10.0, -f64::from(70_000_u32) / 10.0) as f32;
    assert_eq!(far.to_bits(), likelihood_of(u16::MAX).to_bits());
    let record = SiteRecord::new(
        10,
        &[0.5],
        &[
            biallelic(&[0, 70_000, 455], &[1.0, 0.0, 0.0]),
            biallelic(&[30, 0, 300], &[0.0, 1.0, 0.0]),
        ],
    );
    let read: Vec<u32> = record
        .likelihoods(0)
        .map(|l| (l as f32).to_bits())
        .collect();
    assert_eq!(
        read,
        vec![
            1.0_f32.to_bits(),
            0.0_f32.to_bits(),
            likelihood_of(455).to_bits()
        ]
    );
}
