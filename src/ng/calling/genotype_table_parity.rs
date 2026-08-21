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
//! **The oracle is production's own artefact, built by production's own code.**
//! `shape_for(ploidy, n_alleles)` returns the `GenotypeShape` the shipping posterior
//! engine uses, and every field compared below is the one the engine reads. Nothing is
//! re-derived here and nothing is transcribed, so a change to production's enumeration,
//! its fold, or any of its three formulas makes this test fail rather than pass quietly.
//!
//! **Reaching it cost one edit to the frozen tree, and it is the only one.**
//! `posterior_engine.rs` declared `mod shape;` privately, so `GenotypeShape` and
//! `shape_for` could not be named from `src/ng/` whatever their own `pub(crate)`
//! visibility — a call from here failed with ``error[E0603]: module `shape` is
//! private``. The declaration is now `pub(crate)` (owner, 2026-08-21). That is the
//! whole change: no production behaviour moved, nothing was re-exported, and the
//! dependency runs one way — this test reads production, production still names nothing
//! in ng.
//!
//! **This is ng's test; production is only the yardstick**, exactly as
//! [`crate::ng::scanner_parity`] reads `src/ssr/` without writing to it. It is its own
//! file rather than a block inside `genotype_table.rs`'s `mod tests` so that ng's only
//! `use crate::var_calling::` sits in one greppable place and the table's own tests stay
//! free of production.
//!
//! **The oracle stops being production above 255 alleles.** `genotype_order`, which
//! `GenotypeShape::build` enumerates with, iterates `min_allele..(n_alleles as u8)`, so
//! 256 alleles yields no genotypes at all, where the port reaches 65,536. A grid that
//! went that wide would fail on the count assertion and read as though the port were
//! wrong. The widest shape here is 18 alleles.

use std::sync::Arc;

use crate::ng::calling::genotype_table::GenotypeTable;
use crate::ng::types::{AlleleId, Ploidy};
use crate::var_calling::posterior_engine::shape::shape_for;

/// Compare one shape's four quantities against production's `GenotypeShape` and return
/// how many genotypes the **table** holds, so the caller can assert a grid's total
/// reach. The table's own count rather than production's: the two are asserted equal one
/// line earlier, and returning the table's makes the totals evidence about the subject.
///
/// Every comparison is exact. The coefficients are compared by bit pattern rather than
/// with a tolerance, deliberately: the port keeps production's summation order on
/// purpose, so anything but bit-equality means it drifted. Reversing that order alone
/// moves a coefficient by four units in the last place, which any ordinary tolerance
/// would accept.
fn compare_against_production(copies: u8, allele_count: usize) -> usize {
    let ploidy = Ploidy::try_new(copies).expect("the grids start at ploidy 1");
    let table = GenotypeTable::build(ploidy, allele_count);
    let production = shape_for(copies, allele_count);

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
    assert_eq!(
        table.genotype_count(),
        production.n_genotypes,
        "{shape}: count"
    );

    // 2. Every row's allele counts, in production's enumeration order. Compared as one
    //    slice, so a reordering fails as loudly as a wrong count would.
    assert_eq!(
        table.genotype_allele_counts(),
        production.genotype_allele_counts.as_slice(),
        "{shape}: allele counts, or the order they are in"
    );

    // 3. Every log coefficient, to floating-point equality. Compared as one slice of bit
    //    patterns, so a table with the wrong number of rows fails here rather than
    //    having its extra rows go unread.
    let ours: Vec<u64> = table
        .log_multinomial_coeffs()
        .iter()
        .map(|coeff| coeff.to_bits())
        .collect();
    let theirs: Vec<u64> = production
        .log_multinomial_coeffs
        .iter()
        .map(|coeff| coeff.to_bits())
        .collect();
    assert_eq!(
        ours,
        theirs,
        "{shape}: log multinomial coefficients as bit patterns, or how many of them — \
         ours {:?}, production {:?}",
        table.log_multinomial_coeffs(),
        production.log_multinomial_coeffs
    );

    // 4. The homozygous lookup. ng names the allele with an `AlleleId`, production with
    //    a bare `u8`; the widening is the only difference between the two tables.
    let theirs: Vec<Option<AlleleId>> = production
        .homozygous_allele_for
        .iter()
        .map(|entry| entry.map(|allele| AlleleId(u16::from(allele))))
        .collect();
    assert_eq!(
        table.homozygous_alleles(),
        theirs.as_slice(),
        "{shape}: homozygous lookup, or how many entries it has"
    );

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
/// debug profile `cargo test` uses — to compare the same quantities at wider tables.
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
