use super::*;

fn fragment(fragment: u32, allele: u16, error_phred: u8) -> FragmentAllele {
    FragmentAllele {
        sample: 0,
        fragment,
        allele,
        error_phred,
    }
}

fn het(lower: u16, upper: u16) -> Heterozygote {
    Heterozygote::of(lower, upper).expect("two alleles")
}

/// **The thresholds are the logs they are named for**: odds of 199 to 1 and 9 to 1.
#[test]
fn the_thresholds_are_ln_199_and_ln_9() {
    assert_eq!(LINK_MIN, crate::float::ln(199.0));
    assert_eq!(BREAK_MIN, crate::float::ln(9.0));
}

/// **An error is floored at 0.01 and capped at one half** (§2.5): Phred 20 and 40 both read as
/// 0.01, Phred 13 as about 0.05, and Phred 0, 2 and 3 — an error of 1, 0.63 and 0.50 — as one half.
#[test]
fn an_error_is_floored_and_capped() {
    assert_eq!(error_from_phred(40), ERROR_FLOOR);
    assert!((error_from_phred(20) - 0.01).abs() < 1e-15);
    assert!((error_from_phred(13) - 0.050_118_723).abs() < 1e-9);
    for phred in [0, 2, 3] {
        assert_eq!(error_from_phred(phred), 0.5, "Phred {phred}");
    }
}

/// **What one fragment is worth** (§2.5's own arithmetic): at most 3.90 nats however good its
/// bases, below [`LINK_MIN`]; a Phred-13 read allele beside a good one gives 2.77, above
/// [`BREAK_MIN`], so one fragment alone can break a population link.
#[test]
fn one_fragment_is_worth_at_most_3_90_nats() {
    let best = weight(60, 60);
    assert!((best - crate::float::ln(0.9802 / 0.0198)).abs() < 1e-12);
    assert!((best - 3.902).abs() < 1e-3);
    assert!(best < LINK_MIN);
    assert!(2.0 * best >= LINK_MIN, "two good fragments make a link");
    let weak = weight(13, 60);
    assert!((weak - 2.767).abs() < 1e-3, "{weak}");
    assert!(weak >= BREAK_MIN);
}

/// **A read allele at an error of one half weighs exactly nothing, and one at Phred 2 does not
/// reverse a link** (§2.5): without the ceiling, Phred 2's error of 0.63 would vote for the
/// opposite relation.
#[test]
fn a_read_allele_at_one_half_weighs_exactly_nothing() {
    for phred in [0, 2, 3] {
        for other in [0, 13, 30, 60] {
            assert_eq!(weight(phred, other), 0.0, "Phred {phred} beside {other}");
            assert_eq!(weight(other, phred), 0.0);
        }
    }
    let link = score_read_link(
        &[fragment(1, 0, 2), fragment(2, 0, 30), fragment(3, 0, 30)],
        het(0, 1),
        &[fragment(1, 1, 30), fragment(2, 0, 30), fragment(3, 0, 30)],
        het(0, 1),
    );
    assert_eq!(link.log_odds, 2.0 * weight(30, 30));
    assert_eq!(link.relation(), Some(Relation::Cis));
}

/// **Two good fragments that agree link the records; one does not; fragments that disagree
/// cancel; and a fragment showing an allele outside either genotype, or at only one record,
/// says nothing.**
#[test]
fn the_sum_runs_over_fragments_at_both_records_with_alleles_of_both_genotypes() {
    let one = het(0, 1);
    let other = het(0, 2);
    // Fragments 4 and 9 say cis (lower with lower), fragment 6 trans; fragment 7 shows allele 1
    // at the second record, which its genotype 0/2 lacks; fragment 8 is at the first only.
    let at_one = [
        fragment(4, 0, 30),
        fragment(6, 1, 30),
        fragment(7, 0, 30),
        fragment(8, 1, 30),
        fragment(9, 1, 30),
    ];
    let at_other = [
        fragment(4, 0, 30),
        fragment(6, 0, 30),
        fragment(7, 1, 30),
        fragment(9, 2, 30),
    ];
    let w = weight(30, 30);

    let two_and_one = score_read_link(&at_one, one, &at_other, other);
    assert_eq!(two_and_one.fragments, 3);
    assert_eq!(two_and_one.log_odds, w - w + w);
    assert_eq!(
        two_and_one.relation(),
        None,
        "one fragment's worth is not a link"
    );

    let two = score_read_link(&at_one[..1], one, &at_other[..1], other);
    assert_eq!(two.relation(), None);
    let agreeing = score_read_link(
        &[at_one[0], at_one[4]],
        one,
        &[at_other[0], at_other[3]],
        other,
    );
    assert_eq!(agreeing.log_odds, w + w);
    assert_eq!(agreeing.relation(), Some(Relation::Cis));

    let trans = score_read_link(
        &[fragment(1, 0, 30), fragment(2, 1, 30)],
        one,
        &[fragment(1, 2, 30), fragment(2, 0, 30)],
        other,
    );
    assert_eq!(trans.log_odds, -(w + w));
    assert_eq!(trans.relation(), Some(Relation::Trans));
}

/// **A contradiction is judged against the relation it contradicts**: −2.20 nats contradicts
/// cis, not trans, and a link exactly at the threshold counts.
#[test]
fn a_contradiction_is_at_or_beyond_the_break_threshold_against_the_relation() {
    let against_cis = ReadLink {
        log_odds: -BREAK_MIN,
        fragments: 1,
    };
    assert!(against_cis.contradicts(Relation::Cis));
    assert!(!against_cis.contradicts(Relation::Trans));
    let short_of_it = ReadLink {
        log_odds: -BREAK_MIN + 1e-9,
        fragments: 1,
    };
    assert!(!short_of_it.contradicts(Relation::Cis));
    let at_link = ReadLink {
        log_odds: LINK_MIN,
        fragments: 2,
    };
    assert_eq!(at_link.relation(), Some(Relation::Cis));
}

/// **The two orientation rules are each other's inverse**: two records written the same way round
/// are cis, opposite ways trans, and carrying an orientation through a relation gives back one
/// that relation holds with.
#[test]
fn relation_between_and_carry_are_inverses() {
    for one in [false, true] {
        for other in [false, true] {
            let relation = Relation::between(one, other);
            assert_eq!(relation == Relation::Cis, one == other);
            assert_eq!(relation.carry(one), other);
        }
    }
}

/// **No stored Phred is finer than the error floor**: the largest one stored, 60, is an error far
/// below 0.01, so capping Phreds there changes no weight (§2.2). Keeps the two constants coupled.
#[test]
fn the_largest_stored_phred_is_below_the_error_floor() {
    assert!(crate::float::powf(10.0, -f64::from(LARGEST_ERROR_PHRED) / 10.0) < ERROR_FLOOR);
    assert_eq!(weight(LARGEST_ERROR_PHRED, 20), weight(u8::MAX, 20));
}
