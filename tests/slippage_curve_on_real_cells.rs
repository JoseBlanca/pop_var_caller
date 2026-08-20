//! The slippage curve, fitted against the cells two real cohorts actually produced.
//!
//! **The unit tests in `slippage_curve.rs` draw cells on a curve and check the fit finds it
//! back. This one has no known answer to find** — it runs the fit over the per-cell tables the
//! walk wrote for the tomato cohort and for HG002, and asserts the shape numbers reported in
//! `doc/devel/ng/reports/str_slippage_shape_2026-08-20.md` §4.1. It is the step that would catch
//! a fit that is self-consistent and wrong on data.
//!
//! **The two tables were produced with borrowing off**, so every cell speaks from its own tracts
//! alone. Fitting a curve to cells that had already borrowed from their neighbours would be
//! circular, and it would look like a triumph.
//!
//! **They are the ±8 tables, and the ±4 ones they replaced said something different.** The census
//! records a read's length offset over a fixed window and folds the rest into an end bucket; at
//! ±4 that under-measured the level 2.26-fold at 30-repeat homopolymers and made the rise look as
//! though it flattened. Every number asserted below moved when the window widened — see the
//! report's §7.
//!
//! Fixtures: `tests/data/slippage_cells/`. They are the raw `SSR_CELL_TABLE` output of
//! `examples/ng_joint_records_walk.rs`, copied unchanged.

use std::collections::BTreeMap;
use std::path::Path;

use pop_var_caller::ng::parameter_estimation::joint::slippage_curve::{
    FittedCell, RiseShape, SlippageCurveConfig, choose_rise_shape,
};

/// One motif period's contributing cells, read out of a cell table.
///
/// **Only rows the walk marked as fitted are here** — a refused cell has no level, and a cell
/// that borrowed would make the exercise circular. On these two tables nothing borrowed.
fn cells_by_period(path: &Path) -> BTreeMap<u8, Vec<FittedCell>> {
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
    let (period, repeats, reads, fitted, level) = (
        column("period"),
        column("repeats"),
        column("spanning_reads"),
        column("fitted"),
        column("level"),
    );

    let mut by_period: BTreeMap<u8, Vec<FittedCell>> = BTreeMap::new();
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
            .push(FittedCell {
                repeats: field[repeats].parse().expect("a repeat count"),
                level,
                // The count the fitted level says slipped — what sets the cell's own precision.
                slipped_reads: level * spanning,
            });
    }
    for cells in by_period.values_mut() {
        cells.sort_by_key(|cell| cell.repeats);
    }
    by_period
}

fn fixture(name: &str) -> BTreeMap<u8, Vec<FittedCell>> {
    cells_by_period(
        &Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/data/slippage_cells")
            .join(name),
    )
}

fn shape_at(cells: &BTreeMap<u8, Vec<FittedCell>>, period: u8) -> (RiseShape, f64, usize) {
    let curves = choose_rise_shape(
        std::slice::from_ref(&cells[&period]),
        &SlippageCurveConfig::default(),
    )
    .unwrap_or_else(|refusal| panic!("period {period} gave no curve: {refusal:?}"));
    (curves.rise_shape, curves.held_out_error, curves.cells)
}

/// The two cohorts land at **different** shape numbers, which is the finding the fitted family
/// exists for — neither cohort's answer would serve the other.
///
/// **They are no longer opposite, and that correction is the recording window's.** At ±4 tomato
/// fitted 0.00 and HG002 1.00, the two ends of the grid; at ±8 they fit 0.65 and 0.35, both well
/// inside it.
#[test]
fn the_two_cohorts_return_the_shape_numbers_the_report_records() {
    let tomato = fixture("tomato_63_accessions_8mb.csv");
    let hg002 = fixture("hg002_300x_tier.csv");

    let (shape, _, cells) = shape_at(&tomato, 1);
    assert_eq!(cells, 5, "tomato's homopolymer cells");
    assert!(
        (shape.get() - 0.65).abs() < 1e-9,
        "tomato's homopolymers should fit 0.65, got {shape}"
    );

    let (shape, _, cells) = shape_at(&hg002, 1);
    assert_eq!(cells, 23, "HG002's homopolymer cells");
    assert!(
        (shape.get() - 0.35).abs() < 1e-9,
        "HG002's homopolymers should fit 0.35, got {shape}"
    );

    let (shape, _, cells) = shape_at(&hg002, 2);
    assert_eq!(cells, 20, "HG002's dinucleotide cells");
    assert!(
        (shape.get() - 0.70).abs() < 1e-9,
        "HG002's dinucleotides should fit 0.70, got {shape}"
    );
}

/// How far each period's curve lands from a cell it never saw — the number spec §7 reads as the
/// curve's own precision when it weighs a curve against a cell's own answer.
///
/// **These are asserted rather than left to the prose**, so the report's table cannot go stale
/// silently. The bounds are ±0.2 percentage points around what this fit returns, which is the
/// precision the report prints.
#[test]
fn each_periods_curve_lands_where_the_report_says_it_does() {
    let expected = [
        ("tomato_63_accessions_8mb.csv", 1_u8, 5_usize, 0.0734),
        ("hg002_300x_tier.csv", 1, 23, 0.0774),
        ("hg002_300x_tier.csv", 2, 20, 0.1139),
        ("hg002_300x_tier.csv", 3, 4, 0.3199),
        ("hg002_300x_tier.csv", 4, 7, 0.2133),
    ];
    for (name, period, want_cells, want_error) in expected {
        let (_, held_out_error, cells) = shape_at(&fixture(name), period);
        assert_eq!(cells, want_cells, "{name} period {period}");
        assert!(
            (held_out_error - want_error).abs() < 0.002,
            "{name} period {period} predicts a held-out cell to {:.2}%, expected {:.2}%",
            held_out_error * 100.0,
            want_error * 100.0
        );
    }
}

/// **The thin periods are the reason spec §11 doubts the four-cell floor.** HG002's
/// trinucleotides clear it with exactly four cells and predict a held-out cell to 32.0%, nearly
/// three times worse than its dinucleotides do with twenty at 11.4%.
#[test]
fn a_period_at_the_cell_floor_is_far_less_certain_than_a_rich_one() {
    let hg002 = fixture("hg002_300x_tier.csv");
    let (_, at_the_floor, cells) = shape_at(&hg002, 3);
    assert_eq!(cells, SlippageCurveConfig::default().min_cells_for_a_curve);
    let (_, rich, _) = shape_at(&hg002, 2);
    assert!(
        at_the_floor > rich * 2.5,
        "four cells gave {:.1}% and twenty gave {:.1}%",
        at_the_floor * 100.0,
        rich * 100.0
    );
}

/// **The claim that decides the family.** A curve fitted only over 8 to 12 repeats, then asked
/// about 30, must not return a number that is not a probability — which is what an unbounded
/// exponential does here, reaching 205.
#[test]
fn a_curve_fitted_on_a_narrow_window_does_not_explode_outside_it() {
    let hg002 = fixture("hg002_300x_tier.csv");
    let narrow: Vec<FittedCell> = hg002[&1]
        .iter()
        .filter(|cell| (8..=12).contains(&cell.repeats))
        .copied()
        .collect();
    assert_eq!(narrow.len(), 5);

    let curves = choose_rise_shape(&[narrow], &SlippageCurveConfig::default())
        .expect("five cells over four repeat counts");
    let curve = curves.by_group[0].expect("the only group has a line");
    assert_eq!((curve.fitted_from, curve.fitted_to), (8, 12));

    let far = curve.level_at(30);
    assert!(
        (0.0..=1.0).contains(&far),
        "the level at 30 repeats came back as {far}"
    );
    // Held at the 12-repeat end, so it equals what the curve says there.
    assert_eq!(far, curve.level_at(12));
}

/// Every period that clears the floor gives a curve that rises, on both cohorts. A period whose
/// cells will not support one says so rather than returning a falling level.
#[test]
fn every_period_that_clears_the_floor_gives_a_rising_curve() {
    for name in ["tomato_63_accessions_8mb.csv", "hg002_300x_tier.csv"] {
        let cells = fixture(name);
        for (period, at_period) in &cells {
            let outcome = choose_rise_shape(
                std::slice::from_ref(at_period),
                &SlippageCurveConfig::default(),
            );
            match outcome {
                Ok(curves) => {
                    for curve in curves.by_group.iter().flatten() {
                        assert!(
                            curve.slope > 0.0,
                            "{name} period {period} returned a falling curve"
                        );
                        let mut previous = 0.0;
                        for repeats in 1..=60 {
                            let level = curve.level_at(repeats);
                            assert!(
                                level >= previous - 1e-15,
                                "{name} period {period} fell at {repeats} repeats"
                            );
                            assert!((0.0..=1.0).contains(&level));
                            previous = level;
                        }
                    }
                }
                Err(refusal) => assert!(
                    at_period.len() < SlippageCurveConfig::default().min_cells_for_a_curve,
                    "{name} period {period} has {} cells and was refused: {refusal:?}",
                    at_period.len()
                ),
            }
        }
    }
}
