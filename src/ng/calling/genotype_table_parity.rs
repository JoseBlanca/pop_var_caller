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
//! **The oracle is production's own artefact, built by production's own code — frozen since
//! promotion step C5.** Until then this file called `shape_for(ploidy, n_alleles)`, which
//! returned the `GenotypeShape` production's posterior engine used, and compared every field
//! the engine read. Production is being deleted, so each of the 76 shapes the tests compare was
//! written once, at commit `d9e7b076`, to
//! [`testdata/genotype_tables_production.txt`](testdata/genotype_tables_production.txt) — 27,384
//! genotype rows — and the tests compare against that file. Nothing is re-derived here and
//! nothing is transcribed by hand.
//!
//! **This is ng's test; production is only the yardstick.** It is its own file rather than a
//! block inside `genotype_table.rs`'s `mod tests` so that the table's own tests stay free of
//! the comparison.
//!
//! **The oracle stopped being production above 255 alleles.** `genotype_order`, which
//! `GenotypeShape::build` enumerated with, iterated `min_allele..(n_alleles as u8)`, so
//! 256 alleles yielded no genotypes at all, where the port reaches 65,536. The widest shape
//! here is 18 alleles.

use std::collections::HashMap;
use std::sync::{Arc, OnceLock};

use crate::ng::calling::genotype_table::GenotypeTable;
use crate::ng::types::{AlleleId, Ploidy};

/// Production's table for one shape, as read back from the fixture.
struct ProductionShape {
    genotype_count: usize,
    /// Row-major `genotype_count × allele_count`, as production's `genotype_allele_counts`.
    genotype_allele_counts: Vec<u32>,
    /// Bit patterns, as production's `log_multinomial_coeffs`.
    log_multinomial_coeff_bits: Vec<u64>,
    homozygous_allele_for: Vec<Option<AlleleId>>,
}

/// Every shape in the fixture, keyed by `(ploidy, allele count)`, parsed once.
fn production_shapes() -> &'static HashMap<(u8, usize), ProductionShape> {
    static SHAPES: OnceLock<HashMap<(u8, usize), ProductionShape>> = OnceLock::new();
    SHAPES.get_or_init(|| {
        let digit = |c: char| c.to_digit(36).expect("a base-36 allele digit") as usize;
        let mut shapes = HashMap::new();
        let mut lines = include_str!("testdata/genotype_tables_production.txt")
            .lines()
            .filter(|line| !line.starts_with('#'));
        while let Some(header) = lines.next() {
            let fields: Vec<&str> = header.split(' ').collect();
            let ["shape", ploidy, allele_count, genotype_count] = fields[..] else {
                panic!("a shape header, not {header:?}");
            };
            let (ploidy, allele_count, genotype_count): (u8, usize, usize) = (
                ploidy.parse().expect("a ploidy"),
                allele_count.parse().expect("an allele count"),
                genotype_count.parse().expect("a genotype count"),
            );
            let mut shape = ProductionShape {
                genotype_count,
                genotype_allele_counts: vec![0; genotype_count * allele_count],
                log_multinomial_coeff_bits: Vec::with_capacity(genotype_count),
                homozygous_allele_for: Vec::with_capacity(genotype_count),
            };
            for row in 0..genotype_count {
                let line = lines
                    .next()
                    .expect("a genotype line for every genotype counted");
                let fields: Vec<&str> = line.split(' ').collect();
                let [alleles, coeff_bits, homozygous] = fields[..] else {
                    panic!("a genotype line, not {line:?}");
                };
                // Checked here rather than left to the comparison: a digit past the width would
                // land in the next row's counts and fail later under a misleading message.
                assert_eq!(
                    alleles.len(),
                    usize::from(ploidy),
                    "shape {ploidy} {allele_count}: {line:?} does not carry one digit per copy"
                );
                assert!(
                    alleles.as_bytes().is_sorted(),
                    "shape {ploidy} {allele_count}: {line:?} lists its alleles out of order"
                );
                for allele in alleles.chars() {
                    let allele = digit(allele);
                    assert!(
                        allele < allele_count,
                        "shape {ploidy} {allele_count}: {line:?} names an allele past the width"
                    );
                    shape.genotype_allele_counts[row * allele_count + allele] += 1;
                }
                shape
                    .log_multinomial_coeff_bits
                    .push(u64::from_str_radix(coeff_bits, 16).expect("a bit pattern"));
                shape.homozygous_allele_for.push(match homozygous {
                    "-" => None,
                    allele => {
                        assert_eq!(
                            allele.len(),
                            1,
                            "a homozygous allele is one digit: {line:?}"
                        );
                        let allele = digit(allele.chars().next().expect("one digit"));
                        Some(AlleleId(
                            u16::try_from(allele).expect("an allele fits a u16"),
                        ))
                    }
                });
            }
            let previous = shapes.insert((ploidy, allele_count), shape);
            assert!(
                previous.is_none(),
                "shape {ploidy} {allele_count} is frozen twice"
            );
        }
        shapes
    })
}

/// Compare one shape's four quantities against production's frozen `GenotypeShape` and return
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
    let shape = format!("ploidy {copies} over {allele_count} alleles");
    let production = production_shapes()
        .get(&(copies, allele_count))
        .unwrap_or_else(|| panic!("{shape}: production's tables were not frozen for this shape"));

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
        production.genotype_count,
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
    assert_eq!(
        ours,
        production.log_multinomial_coeff_bits,
        "{shape}: log multinomial coefficients as bit patterns, or how many of them — \
         ours {:?}",
        table.log_multinomial_coeffs(),
    );

    // 4. The homozygous lookup. Production named the allele with a bare `u8`; the fixture
    //    widens it to ng's `AlleleId` on reading.
    assert_eq!(
        table.homozygous_alleles(),
        production.homozygous_allele_for.as_slice(),
        "{shape}: homozygous lookup, or how many entries it has"
    );

    table.genotype_count()
}

// ---------------------------------------------------------------------
// The grids
// ---------------------------------------------------------------------

/// The grid the plan names: ploidy 2 and 4, allele counts 1 to 6 — a diploid and a
/// tetraploid locus at every candidate width up to the cap the calling loop will ship,
/// [`DEFAULT_MAX_CANDIDATE_ALLELES`](crate::ng::calling::allele_candidates::DEFAULT_MAX_CANDIDATE_ALLELES),
/// whose value is inherited from production's `DEFAULT_MAX_ALLELES_PER_RECORD` and recorded in
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
