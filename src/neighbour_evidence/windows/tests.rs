use super::*;

/// Records every `spacing` bases from position 1 to `length`.
fn evenly(length: u64, spacing: u64) -> Vec<u64> {
    (0..length / spacing)
        .map(|record| 1 + record * spacing)
        .collect()
}

/// The windows of a whole contig, planned as a streaming caller would: records handed over a few
/// at a time, a plan asked for after each batch.
fn plan_contig(positions: &[u64], samples: usize, sizes: &WindowSizes) -> Vec<WindowPlan> {
    let mut plans = Vec::new();
    let mut read = 0;
    let mut core_start = 0;
    while core_start < positions.len() {
        if read <= core_start {
            read = (read + 997).min(positions.len());
            continue;
        }
        let complete = read == positions.len();
        match plan_window(&positions[..read], core_start, complete, samples, sizes) {
            Some(plan) => {
                core_start = plan.core_end;
                plans.push(plan);
            }
            None => {
                assert!(!complete, "a complete contig always plans");
                read = (read + 997).min(positions.len());
            }
        }
    }
    plans
}

/// **Every record is in exactly one core, in arrival order**, at sizes from full to past the
/// floor (§8.7's invariant at window boundaries).
#[test]
fn cores_cover_every_record_once_and_in_order() {
    let positions = evenly(12_000_000, 40);
    for samples in [1, 57, 1_000, 3_000, 10_000] {
        let plans = plan_contig(&positions, samples, &WindowSizes::default());
        let mut next = 0;
        for plan in &plans {
            assert_eq!(plan.core_start, next, "{samples} samples");
            assert!(plan.core_end > plan.core_start);
            assert!(plan.context_start <= plan.core_start && plan.core_end <= plan.context_end);
            next = plan.core_end;
        }
        assert_eq!(next, positions.len(), "{samples} samples");
    }
}

/// **A small cohort gets the full window**: a 2-Mb core and 1-Mb margins at 57 samples and
/// tomato1's density of 26 records a kilobase.
#[test]
fn a_small_cohort_gets_the_full_window() {
    let positions = evenly(20_000_000, 38);
    let plan = plan_window(
        &positions,
        positions.len() / 2,
        true,
        57,
        &WindowSizes::default(),
    )
    .expect("the contig is complete");
    assert_eq!((plan.core_bases, plan.margin_bases), (2_000_000, 1_000_000));
    assert_eq!(plan.over_budget, 0);
}

/// **At 3,000 samples the window shrinks to 1.5 Mb within the budget** (§6's arithmetic): the
/// core first, to the margins' length, then both to the 500-kb floor.
#[test]
fn three_thousand_samples_shrink_to_the_floor_within_budget() {
    let positions = evenly(20_000_000, 38);
    let plan = plan_window(
        &positions,
        positions.len() / 2,
        true,
        3_000,
        &WindowSizes::default(),
    )
    .expect("the contig is complete");
    assert_eq!((plan.core_bases, plan.margin_bases), (500_000, 500_000));
    assert_eq!(plan.over_budget, 0);
    assert!(plan.records() as u64 * 3_000 <= 120_000_000);
    // 1,300 samples shrink only the core.
    let plan = plan_window(
        &positions,
        positions.len() / 2,
        true,
        1_300,
        &WindowSizes::default(),
    )
    .expect("the contig is complete");
    assert_eq!(plan.margin_bases, 1_000_000);
    assert!(plan.core_bases < 2_000_000 && plan.core_bases >= 1_000_000);
}

/// **Past the floor the budget is exceeded, and the plan says by how much**: margins shorter than
/// 500 kb would cut off the evidence a walk reaches for.
#[test]
fn past_the_floor_the_budget_is_exceeded_and_said() {
    let positions = evenly(20_000_000, 38);
    let plan = plan_window(
        &positions,
        positions.len() / 2,
        true,
        5_000,
        &WindowSizes::default(),
    )
    .expect("the contig is complete");
    assert_eq!((plan.core_bases, plan.margin_bases), (500_000, 500_000));
    assert_eq!(
        plan.over_budget,
        plan.records() as u64 * 5_000 - 120_000_000
    );
}

/// **No plan until the records the window reaches have been read**, and one as soon as they
/// have — or as soon as the contig is known to end.
#[test]
fn a_window_waits_for_its_right_margin() {
    let positions = evenly(10_000_000, 1_000);
    let sizes = WindowSizes::default();
    // 2.5 Mb read: the full core and margin reach 3 Mb.
    assert_eq!(plan_window(&positions[..2_500], 0, false, 10, &sizes), None);
    let plan = plan_window(&positions[..3_001], 0, false, 10, &sizes).expect("3 Mb read");
    assert_eq!((plan.core_end, plan.context_end), (2_000, 3_000));
    let plan = plan_window(&positions[..2_500], 0, true, 10, &sizes).expect("the contig ended");
    assert_eq!((plan.core_end, plan.context_end), (2_000, 2_500));
}
