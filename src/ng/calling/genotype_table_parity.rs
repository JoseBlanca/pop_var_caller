//! **ng's genotype table against production's, value for value.**
//!
//! [`GenotypeTable`](super::genotype_table::GenotypeTable) is a port of production's
//! `GenotypeShape` (`src/var_calling/posterior_engine/shape.rs`), and the thing a port
//! has to prove is not that it is self-consistent but that it agrees with what it was
//! ported from. It holds four quantities per shape — the genotype count, the allele
//! counts in the VCF `PL` order, the log multinomial coefficients, and the homozygous
//! lookup — and **three of the four would be wrong silently** if the port had slipped:
//! a different enumeration order labels every `PL` entry with the wrong genotype, a
//! coefficient off in the last bits shifts every prior, a homozygous entry in the wrong
//! row fires the inbreeding mixture on the wrong genotype. Only a wrong count is loud
//! (`doc/devel/ng/impl_plan/calling_foundations.md`, step C2;
//! `doc/devel/ng/arch/calling_em_loop.md` §8).
//!
//! **This is ng's test; production is only the yardstick.** Nothing here writes to
//! `src/var_calling/`, exactly as [`crate::ng::scanner_parity`] reads `src/ssr/` without
//! writing to it. It is its own file rather than a block inside `genotype_table.rs`'s
//! `mod tests` so that the one `use crate::var_calling::` in all of `src/ng/` sits in
//! one greppable place and the table's own tests stay free of production.
//!
//! ## What the oracle is, and why it is not `GenotypeShape` itself
//!
//! The plan says this test calls `GenotypeShape` directly. **It cannot.**
//! `src/var_calling/posterior_engine.rs` declares `mod shape;` — private — so neither
//! `GenotypeShape` nor `shape_for` can be named from `src/ng/`, whatever their own
//! visibility, and there is no re-export route either (`posterior_engine.rs:61` imports
//! them privately). Compiling a call to `shape_for` from here fails with
//! ``error[E0603]: module `shape` is private``, measured on this branch. The alternative
//! would be widening that declaration, which would edit a frozen tree
//! (`src/ng/mod.rs`: ng "does not edit … `src/var_calling/`") for a test's convenience.
//!
//! So the oracle is assembled here from two things:
//!
//! - **`genotype_order`** — `pub fn` in `pub mod per_group_merger`, reachable, and the
//!   *sole* source of production's enumeration: `collect_non_decreasing` is production's
//!   only enumeration recursion, every consumer reaches it through `genotype_order`, and
//!   this is the one thing here that runs production's own code.
//! - **four transcriptions** — `GenotypeShape::build`'s fold of that enumeration into
//!   the flat count table, plus `log_factorial`, `log_multinomial_coefficient` and
//!   `homozygous_allele`, copied from `shape.rs` including its summation order, since
//!   the coefficients are compared with `==` and not with a tolerance.
//!
//! **What this gives up against calling `GenotypeShape`:** if production ever changed
//! its fold or one of the three formulas, this test would keep passing while ng's table
//! stopped matching production's. All four are closed forms rather than behaviour, which
//! is why the substitution is acceptable — and the **order**, the one part that is a
//! recursion with a comparator rather than a formula, is exactly the part compared
//! against production's own function.
//!
//! **The oracle stops being production above 255 alleles.** `genotype_order` iterates
//! `min_allele..(n_alleles as u8)`, so 256 alleles yields no genotypes at all, where the
//! port reaches 65,536. A grid that went that wide would fail on the count assertion and
//! read as though the port were wrong. The widest shape here is 18 alleles.

use std::sync::Arc;

use crate::ng::calling::genotype_table::GenotypeTable;
use crate::ng::types::{AlleleId, Ploidy};
use crate::var_calling::per_group_merger::genotype_order;

// ---------------------------------------------------------------------
// The oracle: production's `shape.rs`, transcribed
// ---------------------------------------------------------------------

/// Production's fold of `genotype_order` into the flat row-major allele-count table
/// (`src/var_calling/posterior_engine/shape.rs`, `build`).
fn production_allele_counts(copies: u8, n_alleles: usize) -> Vec<u32> {
    let genotypes = genotype_order(copies, n_alleles);
    let n_genotypes = genotypes.len();
    let mut genotype_allele_counts = vec![0_u32; n_genotypes * n_alleles];
    for (g_idx, g) in genotypes.iter().enumerate() {
        let row = &mut genotype_allele_counts[g_idx * n_alleles..(g_idx + 1) * n_alleles];
        for &a in g {
            row[a as usize] += 1;
        }
    }
    genotype_allele_counts
}

/// Production's `log_factorial`, summation order included — the two are compared with
/// `==`, so a different order shows up as a last-bit difference.
///
/// The `as` cast is production's own and is kept for transcription fidelity, not for
/// arithmetic: every `u32` is exactly representable in an `f64`, so `f64::from(i)` — the
/// form the port itself uses — gives the same bits.
fn production_log_factorial(n: u32) -> f64 {
    let mut acc = 0.0_f64;
    for i in 2..=n {
        acc += (i as f64).ln();
    }
    acc
}

/// Production's `log_multinomial_coefficient`. `u32::from(copies)` where production
/// writes `ploidy as u32`; both are the zero-extension of the same eight bits.
fn production_log_multinomial_coefficient(copies: u8, counts: &[u32]) -> f64 {
    let mut log_coef = production_log_factorial(u32::from(copies));
    for &k in counts {
        log_coef -= production_log_factorial(k);
    }
    log_coef
}

/// Production's `homozygous_allele`. Production narrows the result to `u8` on the way
/// into its own table; this keeps `usize` and the caller widens to [`AlleleId`], which
/// no input can distinguish, since `genotype_order` cannot produce an allele index above
/// 254.
fn production_homozygous_allele(counts: &[u32]) -> Option<usize> {
    let mut found: Option<usize> = None;
    for (a, &k) in counts.iter().enumerate() {
        if k == 0 {
            continue;
        }
        match found {
            None => found = Some(a),
            Some(_) => return None,
        }
    }
    found
}

// ---------------------------------------------------------------------
// The comparison
// ---------------------------------------------------------------------

/// Compare one shape's four quantities against the oracle and return how many genotypes
/// the **table** holds, so the caller can assert a grid's total reach. The table's own
/// count rather than the oracle's: the two have just been asserted equal, and returning
/// the table's makes the totals evidence about the subject rather than about the oracle.
///
/// Every comparison is exact. The coefficients use `==` rather than a tolerance
/// deliberately: the port keeps production's summation order on purpose, so anything but
/// bit-equality means it drifted.
fn compare_against_production(copies: u8, allele_count: usize) -> usize {
    let ploidy = Ploidy::try_new(copies).expect("the grids start at ploidy 1");
    let table = GenotypeTable::build(ploidy, allele_count);

    let shape = format!("ploidy {copies} over {allele_count} alleles");

    // 0. The table agrees about which shape it is. Everything below indexes rows by the
    //    width the table declares, so a table holding the right numbers under the wrong
    //    declared width would slice its own rows apart.
    assert_eq!(table.ploidy(), ploidy, "{shape}: the table's own ploidy");
    assert_eq!(
        table.allele_count(),
        allele_count,
        "{shape}: the table's own width"
    );

    // 1. Genotype count.
    let expected_count = genotype_order(copies, allele_count).len();
    assert_eq!(table.genotype_count(), expected_count, "{shape}: count");

    // 2. Every row's allele counts, in production's enumeration order. Compared as one
    //    slice, so a reordering fails as loudly as a wrong count would.
    let expected_counts = production_allele_counts(copies, allele_count);
    assert_eq!(
        table.genotype_allele_counts(),
        expected_counts.as_slice(),
        "{shape}: allele counts, or the order they are in"
    );

    // The other two tables are read a row at a time below, for a failure message that
    // names the genotype — so their lengths are asserted here, or a table carrying extra
    // rows would have them go unread.
    assert_eq!(
        table.log_multinomial_coeffs().len(),
        expected_count,
        "{shape}: one coefficient per genotype"
    );
    assert_eq!(
        table.homozygous_alleles().len(),
        expected_count,
        "{shape}: one homozygous entry per genotype"
    );

    // 3. Every log coefficient, to floating-point equality.
    for (row, counts) in expected_counts.chunks_exact(allele_count).enumerate() {
        let expected = production_log_multinomial_coefficient(copies, counts);
        let got = table.log_multinomial_coeffs()[row];
        assert_eq!(
            got.to_bits(),
            expected.to_bits(),
            "{shape}: row {row} {counts:?} coefficient — got {got}, production {expected}"
        );
    }

    // 4. The homozygous lookup.
    for (row, counts) in expected_counts.chunks_exact(allele_count).enumerate() {
        let expected = production_homozygous_allele(counts).map(|allele| {
            AlleleId(u16::try_from(allele).expect("the grids stay far below 65,536 alleles"))
        });
        assert_eq!(
            table.homozygous_alleles()[row],
            expected,
            "{shape}: row {row} {counts:?} homozygous lookup"
        );
    }

    table.genotype_count()
}

// ---------------------------------------------------------------------
// The grids
// ---------------------------------------------------------------------

/// The grid the plan names: ploidy 2 and 4, allele counts 1 to 6 — a diploid and a
/// tetraploid locus at every candidate width up to the cap the calling loop will ship,
/// `DEFAULT_MAX_CANDIDATE_ALLELES = 6`. That name is not yet a constant in this tree;
/// it is inherited from production's `DEFAULT_MAX_ALLELES_PER_RECORD`
/// (`src/var_calling/per_group_merger.rs`) and recorded in
/// `doc/devel/ng/arch/calling_em_loop.md` §8.
///
/// Every shape here recurs in the wider grid below, so this test's own contribution is
/// its total: it is kept because the plan names this grid, not because it reaches
/// anything the next test does not.
#[test]
fn the_table_matches_production_over_the_diploid_and_tetraploid_grid() {
    let mut genotypes_compared = 0;
    for copies in [2_u8, 4] {
        for allele_count in 1..=6_usize {
            genotypes_compared += compare_against_production(copies, allele_count);
        }
    }
    assert_eq!(
        genotypes_compared, 308,
        "the twelve shapes of the plan's grid hold 308 genotypes between them"
    );
}

/// Wider than the plan asks for: ploidy 1 to 8 — one copy of the genome up to the
/// deepest shape the cache keeps — over allele counts 1 to 8. It adds the odd ploidies,
/// the haploid case, and widths past the candidate cap.
///
/// **It stops at 8 alleles rather than the cache's 16 for cost.** These sixty-four
/// shapes hold 24,301 genotypes; the full 8 × 16 grid holds 2,042,958, eighty-four times
/// as many, and takes the module's tests from hundredths of a second to seconds in the
/// debug profile `cargo test` uses — to compare the same four closed forms at wider
/// tables.
#[test]
fn the_table_matches_production_from_haploid_to_octoploid_up_to_eight_alleles() {
    let mut genotypes_compared = 0;
    for copies in 1..=8_u8 {
        for allele_count in 1..=8_usize {
            genotypes_compared += compare_against_production(copies, allele_count);
        }
    }
    assert_eq!(
        genotypes_compared, 24_301,
        "the sixty-four shapes hold 24,301 genotypes between them"
    );
}

/// **The shapes past the cache bounds, where `build` takes its uncached branch** —
/// ploidy 9 and 10, and 17 and 18 alleles. The two grids above stop at 8 on both axes,
/// which is exactly where the cache stops, so without this the branch every polyploid or
/// wide locus takes would never be compared with production at all.
///
/// That gap is not theoretical. A `log_factorial` capped at 8 leaves the homozygous rows
/// exact — `ln 9! − ln 9!` is zero however the terms are capped — while understating
/// every heterozygous coefficient at ploidy 9 by exactly `ln 9 = 2.197` nats: row `[8, 1]`
/// becomes 0 instead of 2.197. That is a genotype prior tilted toward homozygotes by a
/// factor of nine at every polyploid locus, and nothing else in the suite sees it.
///
/// Small on purpose: the widest shape here is 1,140 genotypes.
#[test]
fn the_table_matches_production_past_the_cache_bounds() {
    let mut genotypes_compared = 0;
    for copies in [9_u8, 10] {
        for allele_count in 1..=4_usize {
            genotypes_compared += compare_against_production(copies, allele_count);
        }
    }
    for copies in [2_u8, 3] {
        for allele_count in [17_usize, 18] {
            genotypes_compared += compare_against_production(copies, allele_count);
        }
    }
    assert_eq!(
        genotypes_compared, 3_083,
        "the twelve shapes past the cache bounds hold 3,083 genotypes between them"
    );
}

/// The cache's own contract, which the plan asks for beside the value comparison because
/// a table that is right but rebuilt at every locus is a different defect from one that
/// is wrong.
///
/// Taken at the **boundary** rather than in the middle of the range, and in both
/// directions: the bound is two comparisons, and turning either from `>` into `>=` stops
/// the cache one shape early — at ploidy 8, at 16 alleles, and at the corner where both
/// meet — which a shape comfortably inside the bounds cannot see.
#[test]
fn a_shape_at_the_cache_bound_is_shared_and_one_past_it_is_not() {
    for (copies, allele_count) in [(8_u8, 16_usize), (2, 16), (8, 1)] {
        let ploidy = Ploidy::try_new(copies).expect("the bound starts at ploidy 1");
        let first = GenotypeTable::build(ploidy, allele_count);
        let second = GenotypeTable::build(ploidy, allele_count);
        assert!(
            Arc::ptr_eq(&first, &second),
            "ploidy {copies} over {allele_count} alleles is inside the cache bounds"
        );
    }
    let past_the_bound = Ploidy::try_new(9).expect("9 is a ploidy");
    let first = GenotypeTable::build(past_the_bound, 2);
    let second = GenotypeTable::build(past_the_bound, 2);
    assert!(
        !Arc::ptr_eq(&first, &second),
        "ploidy 9 is past the bound, so each call builds its own"
    );
    assert_eq!(
        first, second,
        "past the bound the values are still the same"
    );
}
