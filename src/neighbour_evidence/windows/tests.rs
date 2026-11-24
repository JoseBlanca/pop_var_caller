use super::*;

/// The tests' budget: §6's 120 million sample-records, at 40 bytes a sample a record.
const BYTES_A_SAMPLE_RECORD: u64 = 40;

fn sizes() -> WindowSizes {
    WindowSizes {
        budget_bytes: 120_000_000 * BYTES_A_SAMPLE_RECORD,
        ..WindowSizes::default()
    }
}

/// Running totals for `records` records, each holding `samples` samples at
/// [`BYTES_A_SAMPLE_RECORD`], nothing more for the core.
fn uniform(records: usize, samples: usize) -> (Vec<u64>, Vec<u64>) {
    let each = samples as u64 * BYTES_A_SAMPLE_RECORD;
    (
        (0..=records as u64).map(|i| i * each).collect(),
        vec![0; records + 1],
    )
}

/// [`plan_window`] over every sample of every record at [`BYTES_A_SAMPLE_RECORD`].
fn plan_uniform(
    positions: &[u64],
    core_start: usize,
    complete: bool,
    samples: usize,
    sizes: &WindowSizes,
) -> Option<WindowPlan> {
    let (held, core) = uniform(positions.len(), samples);
    plan_window(
        positions,
        WindowBytes {
            held: &held,
            core: &core,
            skipped: None,
            outside: None,
        },
        core_start,
        complete,
        sizes,
    )
}

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
        match plan_uniform(&positions[..read], core_start, complete, samples, sizes) {
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
        let plans = plan_contig(&positions, samples, &sizes());
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
    let plan = plan_uniform(&positions, positions.len() / 2, true, 57, &sizes())
        .expect("the contig is complete");
    assert_eq!((plan.core_bases, plan.margin_bases), (2_000_000, 1_000_000));
    assert_eq!(plan.over_budget, 0);
}

/// **At 3,000 samples the window shrinks to 1.5 Mb within the budget** (§6's arithmetic): the
/// core first, to the margins' length, then both to the 500-kb floor.
#[test]
fn three_thousand_samples_shrink_to_the_floor_within_budget() {
    let positions = evenly(20_000_000, 38);
    let plan = plan_uniform(&positions, positions.len() / 2, true, 3_000, &sizes())
        .expect("the contig is complete");
    assert_eq!((plan.core_bases, plan.margin_bases), (500_000, 500_000));
    assert_eq!(plan.over_budget, 0);
    assert!(plan.records() as u64 * 3_000 <= 120_000_000);
    // 1,300 samples shrink only the core.
    let plan = plan_uniform(&positions, positions.len() / 2, true, 1_300, &sizes())
        .expect("the contig is complete");
    assert_eq!(plan.margin_bases, 1_000_000);
    assert!(plan.core_bases < 2_000_000 && plan.core_bases >= 1_000_000);
}

/// **Past the floor the budget is exceeded, and the plan says by how much**: margins shorter than
/// 500 kb would cut off the evidence a walk reaches for.
#[test]
fn past_the_floor_the_budget_is_exceeded_and_said() {
    let positions = evenly(20_000_000, 38);
    let plan = plan_uniform(&positions, positions.len() / 2, true, 5_000, &sizes())
        .expect("the contig is complete");
    assert_eq!((plan.core_bases, plan.margin_bases), (500_000, 500_000));
    assert_eq!(
        plan.over_budget,
        (plan.records() as u64 * 5_000 - 120_000_000) * BYTES_A_SAMPLE_RECORD
    );
}

/// **No plan until the records the window reaches have been read**, and one as soon as they
/// have — or as soon as the contig is known to end.
#[test]
fn a_window_waits_for_its_right_margin() {
    let positions = evenly(10_000_000, 1_000);
    let sizes = sizes();
    // 2.5 Mb read: the full core and margin reach 3 Mb.
    assert_eq!(
        plan_uniform(&positions[..2_500], 0, false, 10, &sizes),
        None
    );
    let plan = plan_uniform(&positions[..3_001], 0, false, 10, &sizes).expect("3 Mb read");
    assert_eq!((plan.core_end, plan.context_end), (2_000, 3_000));
    let plan = plan_uniform(&positions[..2_500], 0, true, 10, &sizes).expect("the contig ended");
    assert_eq!((plan.core_end, plan.context_end), (2_000, 2_500));
}

/// **The core pays its own cost, the margins do not**: a record's core cost counts only where the
/// record is in the window's core, so the same budget holds a larger window when the cost sits in
/// the core than it would if every record paid it.
#[test]
fn only_the_core_pays_the_core_cost() {
    let positions = evenly(4_000_000, 1_000);
    let records = positions.len();
    let held: Vec<u64> = (0..=records as u64).map(|i| i * 100).collect();
    let core: Vec<u64> = (0..=records as u64).map(|i| i * 900).collect();
    let sizes = WindowSizes {
        budget_bytes: 1_000_000,
        ..WindowSizes::default()
    };
    let bytes = WindowBytes {
        held: &held,
        core: &core,
        skipped: None,
        outside: None,
    };
    let plan = plan_window(&positions, bytes, 1_500, true, &sizes).expect("the contig is complete");
    let core_records = (plan.core_end - plan.core_start) as u64;
    assert_eq!(plan.bytes, plan.records() as u64 * 100 + core_records * 900);
    assert!(plan.bytes <= sizes.budget_bytes);
    // Every record paying 1,000 bytes would allow only 1,000 records.
    assert!(plan.records() > 1_000, "{plan:?}");
}

/// Thinning at some rate; the planner reads only whether there is one.
fn thinned(sizes: WindowSizes) -> WindowSizes {
    WindowSizes {
        thinning: Some(MarginThinning {
            thousandths_per_kb: std::num::NonZeroU32::new(40_000).expect("not zero"),
        }),
        ..sizes
    }
}

/// **A thinned window pays the skipped records' solve state only near its core** (§6, thinned margins): every
/// record holds 100 bytes, of which 60 are its solve state, skipped at every other record. The
/// core and the 50 kb either side pay in full, the far margins 70 bytes a record on average.
#[test]
fn a_thinned_window_pays_skipped_solve_state_only_near_its_core() {
    let positions = evenly(4_000_000, 1_000);
    let records = positions.len();
    let held: Vec<u64> = (0..=records as u64).map(|i| i * 100).collect();
    let core = vec![0; records + 1];
    let skipped: Vec<u64> = (0..=records as u64).map(|i| (i / 2) * 60).collect();
    let sizes = thinned(WindowSizes {
        budget_bytes: u64::MAX,
        ..WindowSizes::default()
    });
    let bytes = WindowBytes {
        held: &held,
        core: &core,
        skipped: Some(&skipped),
        outside: None,
    };
    let plan = plan_window(&positions, bytes, 1_500, true, &sizes).expect("the contig is complete");
    // A 2-Mb core, 1-Mb margins, near margins of 50 kb.
    assert_eq!(
        (plan.context_start, plan.core_end, plan.context_end),
        (500, 3_500, 4_000)
    );
    assert_eq!((plan.near_start, plan.near_end), (1_450, 3_550));
    let near = (plan.near_end - plan.near_start) as u64;
    let far = plan.records() as u64 - near;
    assert_eq!(plan.bytes, near * 100 + far * 70);
    // Without thinning the same costs plan the whole window as walked, at full cost.
    let unthinned = plan_window(
        &positions,
        bytes,
        1_500,
        true,
        &WindowSizes {
            budget_bytes: u64::MAX,
            ..WindowSizes::default()
        },
    )
    .expect("the contig is complete");
    assert_eq!(
        (unthinned.near_start, unthinned.near_end),
        (unthinned.context_start, unthinned.context_end)
    );
    assert_eq!(unthinned.bytes, unthinned.records() as u64 * 100);
}

/// **Thinning lets a window too large unthinned fit its budget**, at a larger size than unthinned
/// would plan.
#[test]
fn thinning_fits_a_larger_window_in_the_same_budget() {
    let positions = evenly(4_000_000, 1_000);
    let records = positions.len();
    let held: Vec<u64> = (0..=records as u64).map(|i| i * 100).collect();
    let core = vec![0; records + 1];
    let skipped: Vec<u64> = (0..=records as u64).map(|i| i * 90).collect();
    let budget = WindowSizes {
        budget_bytes: 140_000,
        ..WindowSizes::default()
    };
    let bytes = WindowBytes {
        held: &held,
        core: &core,
        skipped: Some(&skipped),
        outside: None,
    };
    let unthinned = plan_window(&positions, bytes, 1_500, true, &budget).expect("complete");
    let thinned = plan_window(&positions, bytes, 1_500, true, &thinned(budget)).expect("complete");
    assert!(unthinned.over_budget > 0, "{unthinned:?}");
    assert_eq!(thinned.over_budget, 0, "{thinned:?}");
    assert!(thinned.records() > unthinned.records());
}

/// **`--refine-margin-density` reads records a kb to three decimals**, 0 as no thinning, and
/// writes them back as given.
#[test]
fn a_margin_density_reads_and_writes_to_three_decimals() {
    let read = |text: &str| MarginThinning::parse(text);
    assert_eq!(read("0"), Ok(None));
    assert_eq!(read("0.000"), Ok(None));
    for (text, thousandths, written) in [
        ("40", 40_000, "40"),
        ("0.5", 500, "0.5"),
        (".25", 250, "0.25"),
        ("2.125", 2_125, "2.125"),
    ] {
        let thinning = read(text).expect("a density").expect("not zero");
        assert_eq!(thinning.thousandths_per_kb.get(), thousandths);
        assert_eq!(thinning.to_string(), written);
    }
    for bad in ["", ".", "1.2345", "-1", "1e3", "2000000", "x"] {
        assert!(read(bad).is_err(), "{bad}");
    }
}

/// **Below the margin floor, the core shrinks alone** (`--refine-core-floor`): a window the margin
/// floor cannot fit is planned with its margins at the floor and a shorter core, and it fits.
#[test]
fn a_core_floor_fits_a_window_the_margin_floor_cannot() {
    let positions = evenly(4_000_000, 100);
    let records = positions.len();
    let (held, core) = uniform(records, 1);
    let bytes = WindowBytes {
        held: &held,
        core: &core,
        skipped: None,
        outside: None,
    };
    // A record a 100 bases at 40 bytes: 500-kb margins and a 500-kb core hold 600 kB.
    let budget = 500_000;
    let margin_floor_only = WindowSizes {
        budget_bytes: budget,
        ..WindowSizes::default()
    };
    let at_floor = plan_window(&positions, bytes, 20_000, true, &margin_floor_only)
        .expect("the contig is complete");
    assert!(at_floor.over_budget > 0, "{at_floor:?}");
    assert_eq!(
        (at_floor.core_bases, at_floor.margin_bases),
        (500_000, 500_000)
    );

    let with_core_floor = WindowSizes {
        core_floor_bases: 100_000,
        ..margin_floor_only
    };
    let plan = plan_window(&positions, bytes, 20_000, true, &with_core_floor)
        .expect("the contig is complete");
    assert_eq!(plan.over_budget, 0, "{plan:?}");
    assert_eq!(plan.margin_bases, 500_000);
    // The longest core in steps of 50 kb that fits: 1 Mb of margins is 400 kB, leaving 100 kB.
    assert_eq!(plan.core_bases, 250_000);
}

/// **The core floor is a floor, not a size**: with the budget to spare, the core is as long as
/// ever, and a core floor at the margin floor plans what no core floor plans.
#[test]
fn a_core_floor_changes_only_windows_over_the_margin_floor() {
    let positions = evenly(4_000_000, 100);
    let (held, core) = uniform(positions.len(), 1);
    let bytes = WindowBytes {
        held: &held,
        core: &core,
        skipped: None,
        outside: None,
    };
    for budget in [u64::MAX, 3_000_000, 1_300_000, 600_000] {
        let base = WindowSizes {
            budget_bytes: budget,
            ..WindowSizes::default()
        };
        let lower = WindowSizes {
            core_floor_bases: 50_000,
            ..base
        };
        assert_eq!(
            plan_window(&positions, bytes, 20_000, true, &base),
            plan_window(&positions, bytes, 20_000, true, &lower),
            "budget {budget}"
        );
    }
}

/// **A window pays for the records the caller holds outside it** (code review 2026-11-24, M1):
/// 100 records a kb over 3 Mb, the core starting at 1.5 Mb, a caller holding every record. With
/// each record's outside cost given, a window's bytes are what it spans plus that cost of every
/// record before and after it; at full size that is enough to push it over a budget it fit
/// without, so it shrinks: when written, to an 800-kb core and margins.
#[test]
fn a_window_pays_for_the_records_held_outside_it() {
    let positions = evenly(3_000_000, 10);
    let core_start = 150_000;
    let (held, core) = uniform(positions.len(), 1);
    let outside: Vec<u64> = (0..=positions.len() as u64).map(|i| i * 10).collect();
    let sizes = WindowSizes {
        // The full window spans records 50,000 to the end, 250,000 of them, and just fits.
        budget_bytes: 250_000 * BYTES_A_SAMPLE_RECORD + 1_000,
        ..WindowSizes::default()
    };
    let plan_with = |outside: Option<&[u64]>| {
        plan_window(
            &positions,
            WindowBytes {
                held: &held,
                core: &core,
                skipped: None,
                outside,
            },
            core_start,
            true,
            &sizes,
        )
        .expect("the contig is complete")
    };
    let without = plan_with(None);
    assert_eq!(
        (without.core_bases, without.margin_bases),
        (2_000_000, 1_000_000)
    );
    assert_eq!(without.outside_bytes, 0);
    let with = plan_with(Some(&outside));
    assert_eq!(
        (with.core_bases, with.margin_bases),
        (800_000, 800_000),
        "{with:?}"
    );
    let records_outside = (with.context_start + positions.len() - with.context_end) as u64;
    assert_eq!(with.outside_bytes, records_outside * 10);
    assert_eq!(
        with.bytes,
        held[with.context_end] - held[with.context_start] + with.outside_bytes
    );
}
