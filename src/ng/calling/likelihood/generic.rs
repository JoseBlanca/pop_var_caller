//! The SNP/indel closed form — how probable one sample's reads are at an ordinary site,
//! given each candidate genotype.
//!
//! Spec §3 in code. This file starts with the one piece the formula needs before it can be
//! written: **where a wrong read's probability goes.**

use crate::ng::calling::{CandidateAlleles, GenotypeIdx, GenotypeTableView};
use crate::ng::locus_generation::LocusKind;
use crate::ng::types::AlleleId;

/// How many bases a misread could have gone to — **three, and it is a physical fact rather
/// than a tuning choice**: three bases to go wrong into, one to come back to.
///
/// Named because it is kept in sync with three other places: spec §3.5, the parameter
/// pre-pass's own noise model, and whatever the row function does with `log` of it.
pub const ERROR_SPREAD_BASES: f64 = 3.0;

/// The divisor where the model has nothing to say — **one, meaning the error mass is left
/// unspread**, which is the conservative direction and favours the reference.
pub const NO_ERROR_SPREAD: f64 = 1.0;

/// How many things a read could have shown, given that it is wrong — `m(a, g)` in spec §3.5.
///
/// `3.0` where the observation differs from **every** allele the genotype carries by a
/// substitution at exactly one position, `1.0` otherwise.
///
/// # Why a divisor at all
///
/// A read the genotype cannot produce is wrong. **But wrong how?** If the individual carries
/// `A` at this base and the read shows `C`, the chance of *that particular* misread is not the
/// chance of any misread — there were three bases it could have gone to. Dividing the error
/// mass by three is the physical fact, and it is what the parameter pre-pass's own noise model
/// already assumes: *three bases to go wrong into, one to come back to*.
///
/// **The size is `log 3` per wrongly-explained read — 1.10 nats, 4.8 Phred — and it does not
/// cancel**, because how many reads a genotype calls wrong varies by genotype. Dividing by
/// three makes a wrong read *less* probable, so calling a read wrong costs more, and the
/// divisor therefore **favours the heterozygote**.
///
/// The three vendored callers disagree about this and so do the two halves of production:
/// GATK divides by three, freebayes and production's SNP path divide by nothing, and
/// production's own STR substitution term divides by three.
///
/// # Why `1.0` everywhere else, and why *every* allele
///
/// The second case covers insertions, deletions and multi-position differences, where there is
/// **no finite set of things a wrong read could have shown** and any divisor would be invented.
/// Leaving the mass unspread is the conservative choice, and conservative here means favouring
/// the reference — the direction a caller should err in when its model runs out.
///
/// **"Every allele the genotype carries", not "some":** a read one substitution from one of a
/// heterozygote's alleles and an insertion away from the other has no clean three-way spread
/// either, so it takes the conservative divisor. That is what makes `m` a property of the
/// `(allele, genotype)` pair rather than of the allele pair alone.
///
/// # Why it is computed once per locus
///
/// **`m` is a property of the allele pair and not of the read**, so it costs nothing per read:
/// the table is filled once per locus over the projected sequences the merge already unified,
/// and the row function reads it. That unification is by exact byte match, which is sound only
/// because indels were left-aligned upstream — the same reason two samples showing one deletion
/// land on one allele.
///
/// # The layout
///
/// `out` is `genotype_count × allele_count`, **row-major by genotype** — the same shape and the
/// same order as [`GenotypeTableView::genotype_allele_counts`], so a reader holding one row of
/// counts holds the matching row of divisors at the same offset. [`DivisorTable`] is the one
/// way to read it, and it carries the stride so that no caller has to supply one.
///
/// # Panics
///
/// **In release as well as debug**, on three caller bugs (spec §8): an `out` whose length is
/// not exactly `genotype_count × allele_count`; a genotype table whose allele count disagrees
/// with the candidate table's, which is a table built for a different locus and would divide
/// the wrong reads by three with nothing saying so; and a locus that is not
/// [`LocusKind::Generic`], whose substitution term is a different rate on a different model
/// (spec §4.3).
///
/// **The length check is an equality and not "at least enough", which matters twice over.** A
/// longer buffer leaves its tail unwritten, so every genotype past the first reads a slot the
/// fill never touched — a real number and a wrong one. And it is what makes [`DivisorTable`]'s
/// own bound meaningful: a longer table admits a genotype index that should have been out of
/// range.
pub fn fill_error_spread_divisors(
    alleles: &CandidateAlleles,
    genotypes: &GenotypeTableView<'_>,
    out: &mut [f64],
) {
    // **A repeat tract has no business here.** Its substitution term is a different rate on a
    // different model (spec §4.3), and this divisor describes neither — so a table filled for
    // one would be a number with no meaning quietly reaching the wrong row builder.
    assert!(
        matches!(alleles.kind(), LocusKind::Generic),
        "the error-spread divisor is the SNP/indel path's; this locus is {:?}",
        alleles.kind()
    );
    let allele_count = alleles.len();
    assert_eq!(
        allele_count,
        genotypes.allele_count(),
        "the genotype table is built over {} alleles and the candidate table holds {}, so one \
         of them belongs to a different locus",
        genotypes.allele_count(),
        allele_count
    );
    assert_eq!(
        out.len(),
        genotypes.genotype_count() * allele_count,
        "the divisor table needs one entry per (genotype, allele) — {} genotypes × {} alleles \
         = {}, not {}",
        genotypes.genotype_count(),
        allele_count,
        genotypes.genotype_count() * allele_count,
        out.len()
    );

    // One substitution apart is a property of the allele pair, so it is answered once here and
    // looked up below rather than recomputed for every genotype that carries the allele.
    let mut one_substitution_apart = vec![false; allele_count * allele_count];
    for (left, left_bases) in alleles.iter().enumerate() {
        for (right, right_bases) in alleles.iter().enumerate() {
            one_substitution_apart[left * allele_count + right] =
                differ_by_one_substitution(left_bases, right_bases);
        }
    }

    let counts = genotypes.genotype_allele_counts();
    for genotype in 0..genotypes.genotype_count() {
        let carried_copies = &counts[genotype * allele_count..(genotype + 1) * allele_count];
        for observed_allele in 0..allele_count {
            // An allele the genotype carries differs from itself at zero positions, never at
            // exactly one, so it falls out of the `all` below without a special case — and the
            // divisor is never read for an observation the genotype explains anyway.
            let every_carried_is_one_substitution_away = carried_copies
                .iter()
                .enumerate()
                .filter(|&(_, &copies)| copies > 0)
                .all(|(carried_allele, _)| {
                    one_substitution_apart[observed_allele * allele_count + carried_allele]
                });
            out[genotype * allele_count + observed_allele] =
                if every_carried_is_one_substitution_away {
                    ERROR_SPREAD_BASES
                } else {
                    NO_ERROR_SPREAD
                };
        }
    }
}

/// A filled divisor table, with the stride it was filled at.
///
/// **The stride travels with the buffer rather than being handed in at each lookup, and that
/// is the whole reason this type exists.** An accessor taking `(values, allele_count, genotype,
/// allele)` cannot check that `allele_count` is the stride the buffer was actually filled at —
/// so reading a three-allele table at a stride of two returns a real divisor from the wrong
/// row, on half the lookups, with nothing to panic about. Measured on one three-allele diploid
/// locus: six of twelve lookups silently disagree.
///
/// That is exactly the failure this step exists to prevent — `log 3` in the wrong direction and
/// nothing crashes — so the fix is structural rather than another assertion. The crate argues
/// the same case against itself twice already: [`CandidateAlleles::bases_of`] returns an
/// `Option` because indexing "would hand back a real but wrong allele without complaint", and
/// [`GenotypeIdx`] carries the same warning about rows meaning different genotypes at different
/// shapes.
///
/// [`CandidateAlleles::bases_of`]: super::super::CandidateAlleles::bases_of
#[derive(Copy, Clone, Debug)]
pub struct DivisorTable<'a> {
    values: &'a [f64],
    allele_count: usize,
}

impl<'a> DivisorTable<'a> {
    /// Wrap a buffer [`fill_error_spread_divisors`] filled, against the genotype table it was
    /// filled for.
    ///
    /// **The genotype view is the argument rather than a bare stride**, so the two dimensions
    /// come from the same place the fill got them and cannot be supplied separately.
    ///
    /// # Panics
    ///
    /// **In release as well as debug**, on a buffer whose length is not
    /// `genotype_count × allele_count` — the same check the fill makes, repeated here because a
    /// caller can wrap any slice.
    #[must_use]
    pub fn over(values: &'a [f64], genotypes: &GenotypeTableView<'_>) -> Self {
        assert_eq!(
            values.len(),
            genotypes.genotype_count() * genotypes.allele_count(),
            "a divisor table for {} genotypes over {} alleles holds {} entries, not {}",
            genotypes.genotype_count(),
            genotypes.allele_count(),
            genotypes.genotype_count() * genotypes.allele_count(),
            values.len()
        );
        Self {
            values,
            allele_count: genotypes.allele_count(),
        }
    }

    /// One genotype's whole row of divisors — one entry per allele, in allele order.
    ///
    /// **The shape the row function wants**, because its inner loop already holds the matching
    /// row of copy counts as a slice and walks the two together.
    ///
    /// `None` where this table holds no such genotype, for the reason [`GenotypeIdx`] gives:
    /// row 4 of a triallelic diploid table and row 4 of a tetraploid one are different
    /// genotypes, so an index from another shape must not quietly resolve.
    #[must_use]
    pub fn row_of(&self, genotype: GenotypeIdx) -> Option<&'a [f64]> {
        let start = (genotype.get() as usize).checked_mul(self.allele_count)?;
        self.values.get(start..start + self.allele_count)
    }

    /// The divisor for one `(genotype, allele)` pair.
    ///
    /// # Panics
    ///
    /// **In release as well as debug**, on a pair this table does not hold — a genotype from
    /// another shape, or an allele id minted at another locus, which is exactly the case
    /// [`AlleleId`]'s own documentation says is caught when the table is read.
    #[must_use]
    pub fn at(&self, genotype: GenotypeIdx, allele: AlleleId) -> f64 {
        let allele = usize::from(allele.get());
        assert!(
            allele < self.allele_count,
            "allele {allele} is past the {} this locus is called over",
            self.allele_count
        );
        let row = self.row_of(genotype).unwrap_or_else(|| {
            panic!(
                "genotype {} is past the {} this table was filled for",
                genotype.get(),
                self.values.len() / self.allele_count
            )
        });
        row[allele]
    }

    /// How many alleles the locus is called over — the table's stride.
    #[must_use]
    pub fn allele_count(&self) -> usize {
        self.allele_count
    }
}

/// Whether two projected allele sequences differ by a substitution at exactly one position.
///
/// **Different lengths are never a substitution, and the early return is what says so.** An
/// insertion or a deletion changes how many bases there are, and no count of differing
/// positions describes it — which is the whole of why those get the conservative divisor.
/// Without that return the comparison below would `zip`, which **truncates to the shorter
/// sequence**, so `ACGT` against `ATG` would come back as one differing position and be called
/// a substitution. (It would not go wrong on a prefix or a suffix, where truncation leaves zero
/// differences — which is why a fixture built only from those cannot see the check disappear.)
///
/// Two identical sequences differ at zero positions, which is not exactly one, so an allele is
/// never one substitution from itself.
fn differ_by_one_substitution(left: &[u8], right: &[u8]) -> bool {
    if left.len() != right.len() {
        return false;
    }
    let mut differing = 0usize;
    for (a, b) in left.iter().zip(right) {
        if a != b {
            differing += 1;
            if differing > 1 {
                return false;
            }
        }
    }
    differing == 1
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::calling::GenotypeTable;
    use crate::ng::locus_generation::LocusKind;
    use crate::ng::types::Ploidy;

    /// A locus over the bases given, the first of them the reference.
    fn locus(alleles: &[&[u8]]) -> CandidateAlleles {
        let mut table = CandidateAlleles::new(alleles[0].into(), LocusKind::Generic);
        for bases in &alleles[1..] {
            table.admit((*bases).into());
        }
        table
    }

    fn diploid(allele_count: usize) -> std::sync::Arc<GenotypeTable> {
        GenotypeTable::build(Ploidy::try_new(2).expect("two is a ploidy"), allele_count)
    }

    /// Fill the table, so a test reads it through [`DivisorTable`] rather than open-coding the
    /// stride.
    fn divisors(alleles: &CandidateAlleles, table: &GenotypeTable) -> Vec<f64> {
        let view = table.view();
        let mut out = vec![0.0; view.genotype_count() * alleles.len()];
        fill_error_spread_divisors(alleles, &view, &mut out);
        out
    }

    /// The filled buffer, read through the type that carries its stride.
    fn table_over<'a>(out: &'a [f64], table: &GenotypeTable) -> DivisorTable<'a> {
        DivisorTable::over(out, &table.view())
    }

    /// Which genotype index carries exactly these copies, so a test can name a genotype by what
    /// it is rather than by an index the table's ordering happens to give it.
    fn genotype_carrying(table: &GenotypeTable, copies: &[u32]) -> GenotypeIdx {
        let view = table.view();
        let allele_count = view.allele_count();
        (0..view.genotype_count())
            .find(|genotype| {
                &view.genotype_allele_counts()
                    [genotype * allele_count..(genotype + 1) * allele_count]
                    == copies
            })
            .map(|genotype| GenotypeIdx(genotype as u32))
            .expect("the fixture names a genotype the table holds")
    }

    // ---- the one-substitution predicate ----

    #[test]
    fn one_differing_base_is_a_substitution_and_two_are_not() {
        assert!(differ_by_one_substitution(b"ACGT", b"ACCT"));
        assert!(!differ_by_one_substitution(b"ACGT", b"ATCT"));
        assert!(!differ_by_one_substitution(b"ACGT", b"ACGT"));
    }

    /// An insertion or a deletion is not a substitution however few bases it moves, because no
    /// count of differing positions describes it.
    ///
    /// **The last pair is the one that guards the length check**, and the first three cannot.
    /// The comparison `zip`s, and `zip` truncates to the shorter sequence — so on a prefix or a
    /// suffix relation it sees zero differing positions and answers `false` whether the length
    /// check is there or not. `ACGT` against `ATG` truncates to three positions differing at
    /// one, so without the check it would come back a substitution.
    #[test]
    fn an_indel_is_never_one_substitution_however_short() {
        assert!(!differ_by_one_substitution(b"AC", b"ACGT"));
        assert!(!differ_by_one_substitution(b"ACGT", b"ACG"));
        assert!(!differ_by_one_substitution(b"A", b"AC"));
        assert!(!differ_by_one_substitution(b"ACGT", b"ATG"));
        assert!(!differ_by_one_substitution(b"ATG", b"ACGT"));
        assert!(!differ_by_one_substitution(b"AT", b"ACGT"));
        assert!(!differ_by_one_substitution(b"ACG", b"ACCT"));
    }

    /// The predicate does not care which way round it is asked.
    ///
    /// **This carries more weight than it looks.** The pair table is filled over every ordered
    /// pair, so its layout is unobservable *because* the predicate is symmetric — transposing
    /// either the write or the read is a no-op, and nothing else in the file would notice. This
    /// test is what makes that safe rather than lucky.
    ///
    /// **The unequal-length pairs differ inside the overlap on purpose.** A prefix pair like
    /// `AC` against `ACGT` answers `false` both ways under a one-sided length test — one written
    /// as *left longer than right* rather than *lengths differ* — as well as under the real one,
    /// so a fixture built only from prefixes cannot see the asymmetry it is named for.
    #[test]
    fn the_predicate_is_symmetric() {
        for (left, right) in [
            (&b"ACGT"[..], &b"ACCT"[..]),
            (&b"ACGT"[..], &b"ATCT"[..]),
            (&b"AC"[..], &b"ACGT"[..]),
            (&b"AT"[..], &b"ACGT"[..]),
            (&b"ACG"[..], &b"ACCT"[..]),
            (&b"A"[..], &b"CG"[..]),
        ] {
            assert_eq!(
                differ_by_one_substitution(left, right),
                differ_by_one_substitution(right, left),
                "{} against {}",
                String::from_utf8_lossy(left),
                String::from_utf8_lossy(right)
            );
        }
    }

    // ---- the three classes spec §3.5 names ----

    /// **The substitution class.** A biallelic SNP: the reference homozygote cannot explain the
    /// alternative read, and the two alleles differ at exactly one position, so the alternative
    /// gets the three-way spread.
    #[test]
    fn a_read_one_substitution_from_the_only_carried_allele_gets_the_three_way_spread() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let out = divisors(&alleles, &table);
        let hom_ref = genotype_carrying(&table, &[2, 0]);

        assert_eq!(table_over(&out, &table).at(hom_ref, AlleleId(1)), 3.0);
    }

    /// **The multi-position class.** Two alleles the same length differing at two positions have
    /// no three-way spread: there is no single base that went wrong.
    #[test]
    fn a_read_differing_at_two_positions_gets_no_spread() {
        let alleles = locus(&[b"ACGT", b"ATCT"]);
        let table = diploid(2);
        let out = divisors(&alleles, &table);
        let hom_ref = genotype_carrying(&table, &[2, 0]);

        assert_eq!(table_over(&out, &table).at(hom_ref, AlleleId(1)), 1.0);
    }

    /// **The indel class.** A deletion has no finite set of things a wrong read could have
    /// shown, so the mass is left unspread — the conservative choice, which favours the
    /// reference.
    #[test]
    fn a_read_carrying_an_indel_gets_no_spread() {
        let alleles = locus(&[b"ACGT", b"AT"]);
        let table = diploid(2);
        let out = divisors(&alleles, &table);
        let hom_ref = genotype_carrying(&table, &[2, 0]);

        assert_eq!(table_over(&out, &table).at(hom_ref, AlleleId(1)), 1.0);
    }

    /// **Every carried allele, not some** — the case that distinguishes the rule from the
    /// looser one, and the reason `m` is a property of the pair and not of the allele.
    ///
    /// A heterozygote carrying the reference and a deletion, against an observation one
    /// substitution from the reference: one carried allele is a substitution away, the other is
    /// an indel away, so the spread is refused. Under an `any` rule this would be 3.0, and the
    /// two single-allele tests above would pass either way.
    #[test]
    fn a_genotype_carrying_an_indel_refuses_the_spread_for_every_observation() {
        let alleles = locus(&[b"ACGT", b"ACCT", b"AT"]);
        let table = diploid(3);
        let out = divisors(&alleles, &table);
        let het_ref_and_deletion = genotype_carrying(&table, &[1, 0, 1]);
        let hom_ref = genotype_carrying(&table, &[2, 0, 0]);

        // Allele 1 is one substitution from allele 0 and a different length from allele 2.
        assert_eq!(
            table_over(&out, &table).at(het_ref_and_deletion, AlleleId(1)),
            1.0
        );
        // …and against the reference homozygote alone it does get the spread, so the fixture is
        // not simply one where nothing ever would.
        assert_eq!(table_over(&out, &table).at(hom_ref, AlleleId(1)), 3.0);
    }

    /// A heterozygote whose two alleles are *both* one substitution from the observation does
    /// get the spread — the other side of the same rule.
    #[test]
    fn a_genotype_whose_every_allele_is_one_substitution_away_gets_the_spread() {
        // Three alleles at one base: the observation `G` is one substitution from both `A` and
        // `C`.
        let alleles = locus(&[b"A", b"C", b"G"]);
        let table = diploid(3);
        let out = divisors(&alleles, &table);
        let het = genotype_carrying(&table, &[1, 1, 0]);

        assert_eq!(table_over(&out, &table).at(het, AlleleId(2)), 3.0);
    }

    /// An allele the genotype carries never gets the spread, because it differs from itself at
    /// zero positions rather than at one. The divisor is not read for such an observation — it
    /// is on the explained side of the formula — but a table that gave it 3.0 would mean the
    /// predicate had been written as *at most* one rather than *exactly* one, which would also
    /// give an identical pair of alleles the spread.
    #[test]
    fn an_allele_the_genotype_carries_is_not_one_substitution_from_itself() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let out = divisors(&alleles, &table);
        let hom_ref = genotype_carrying(&table, &[2, 0]);
        let het = genotype_carrying(&table, &[1, 1]);

        assert_eq!(table_over(&out, &table).at(hom_ref, AlleleId(0)), 1.0);
        assert_eq!(table_over(&out, &table).at(het, AlleleId(0)), 1.0);
        assert_eq!(table_over(&out, &table).at(het, AlleleId(1)), 1.0);
    }

    // ---- shape and layout ----

    /// The table is genotype-major, matching the genotype table's own counts row for row. A
    /// fill written allele-major produces the same *set* of numbers, so only reading a
    /// specific pair catches the transposition — and the fixture is deliberately not square,
    /// because a square one cannot.
    #[test]
    fn the_table_is_genotype_major_and_matches_the_counts_row_for_row() {
        // Three alleles, so six diploid genotypes: not square, and 3 ≠ 6.
        let alleles = locus(&[b"A", b"C", b"AT"]);
        let table = diploid(3);
        let view = table.view();
        let out = divisors(&alleles, &table);

        assert_eq!(out.len(), 6 * 3);
        for genotype in 0..view.genotype_count() {
            let carried = &view.genotype_allele_counts()[genotype * 3..(genotype + 1) * 3];
            for observed in 0..3 {
                let expected = if carried
                    .iter()
                    .enumerate()
                    .filter(|&(_, &copies)| copies > 0)
                    .all(|(carried_allele, _)| {
                        differ_by_one_substitution(
                            alleles
                                .bases_of(AlleleId(observed as u16))
                                .expect("an allele"),
                            alleles
                                .bases_of(AlleleId(carried_allele as u16))
                                .expect("an allele"),
                        )
                    }) {
                    3.0
                } else {
                    1.0
                };
                assert_eq!(
                    table_over(&out, &table)
                        .at(GenotypeIdx(genotype as u32), AlleleId(observed as u16)),
                    expected,
                    "genotype {genotype}, allele {observed}"
                );
                // …and the row accessor agrees with the pair accessor, which is what lets the
                // row function walk a genotype's divisors beside its copy counts.
                assert_eq!(
                    table_over(&out, &table)
                        .row_of(GenotypeIdx(genotype as u32))
                        .expect("the table holds this genotype")[observed],
                    expected
                );
            }
        }
    }

    /// Ploidy is not in the rule: what a genotype carries is what matters, not how many copies
    /// of it. A tetraploid carrying one copy of an indel refuses the spread exactly as a diploid
    /// heterozygote does.
    #[test]
    fn the_rule_reads_which_alleles_are_carried_and_not_how_many_copies() {
        let alleles = locus(&[b"ACGT", b"ACCT", b"AT"]);
        let tetraploid = GenotypeTable::build(Ploidy::try_new(4).expect("four is a ploidy"), 3);
        let view = tetraploid.view();
        let mut out = vec![0.0; view.genotype_count() * 3];
        fill_error_spread_divisors(&alleles, &view, &mut out);

        let three_ref_one_deletion = genotype_carrying(&tetraploid, &[3, 0, 1]);
        let four_ref = genotype_carrying(&tetraploid, &[4, 0, 0]);

        assert_eq!(
            table_over(&out, &tetraploid).at(three_ref_one_deletion, AlleleId(1)),
            1.0
        );
        assert_eq!(table_over(&out, &tetraploid).at(four_ref, AlleleId(1)), 3.0);
    }

    /// **A locus with only its reference**, which is what `CandidateAlleles::new` produces and
    /// the commonest shape in a genome. One genotype at any ploidy, one entry, and the divisor
    /// is 1.0 — the reference is not one substitution from itself. Nothing else here goes below
    /// two alleles.
    #[test]
    fn a_locus_with_only_its_reference_has_one_entry_and_no_spread() {
        let alleles = locus(&[b"ACGT"]);
        for copies in [1u8, 2, 4] {
            let table = GenotypeTable::build(Ploidy::try_new(copies).expect("a fixture ploidy"), 1);
            let view = table.view();
            assert_eq!(view.genotype_count(), 1);

            let mut out = vec![0.0; 1];
            fill_error_spread_divisors(&alleles, &view, &mut out);

            assert_eq!(
                table_over(&out, &table).at(GenotypeIdx(0), AlleleId::REFERENCE),
                1.0
            );
        }
    }

    /// **At the top of the range**, where nothing else in this file reaches: six alleles, ng's
    /// own default cap, and sixteen, the ceiling production refuses to be configured above. The
    /// point is not the arithmetic — the classes are already pinned — but that the fill's own
    /// shape holds at `21` and `136` genotypes, where an off-by-one in either loop bound or in
    /// the stride has room to show.
    #[test]
    fn the_fill_holds_its_shape_at_the_top_of_the_allele_range() {
        for allele_count in [6usize, 16] {
            // Every allele one substitution from every other: one shared prefix, one varying
            // base, so every pair differs at exactly one position.
            let bases: Vec<Vec<u8>> = (0..allele_count)
                .map(|at| {
                    let mut spelled = b"ACGTACGTAC".to_vec();
                    spelled.push(b'A' + at as u8);
                    spelled
                })
                .collect();
            let refs: Vec<&[u8]> = bases.iter().map(Vec::as_slice).collect();
            let alleles = locus(&refs);
            let table = diploid(allele_count);
            let view = table.view();
            let mut out = vec![f64::NAN; view.genotype_count() * allele_count];

            fill_error_spread_divisors(&alleles, &view, &mut out);

            let divisors_table = table_over(&out, &table);
            assert_eq!(divisors_table.allele_count(), allele_count);
            // A homozygote's own allele is never one substitution from itself; every other
            // allele is one substitution from it, so every other allele gets the spread.
            let hom_ref_counts: Vec<u32> = std::iter::once(2)
                .chain(std::iter::repeat_n(0, allele_count - 1))
                .collect();
            let hom_ref = genotype_carrying(&table, &hom_ref_counts);
            let row = divisors_table
                .row_of(hom_ref)
                .expect("a genotype the table holds");
            assert_eq!(row[0], 1.0);
            assert!(row[1..].iter().all(|&divisor| divisor == 3.0));
        }
    }

    /// **Every cell is written, whatever the buffer held.** The other tests all pass a freshly
    /// zeroed buffer, where a cell the fill skipped reads as `0.0` — neither divisor, so it
    /// would be caught, but by luck rather than on purpose. Poisoning the buffer with `NaN`
    /// makes a skipped cell fail deliberately, and the buffer is caller scratch reused across
    /// loci, so what it held is the previous locus's answer rather than zero.
    #[test]
    fn every_cell_is_written_over_whatever_the_buffer_held() {
        let alleles = locus(&[b"ACGT", b"ACCT", b"AT"]);
        let table = diploid(3);
        let view = table.view();
        let mut out = vec![f64::NAN; view.genotype_count() * 3];

        fill_error_spread_divisors(&alleles, &view, &mut out);

        for (cell, value) in out.iter().enumerate() {
            assert!(
                *value == 1.0 || *value == 3.0,
                "cell {cell} was left holding {value}"
            );
        }
    }

    #[test]
    #[should_panic(expected = "one entry per (genotype, allele)")]
    fn a_buffer_too_short_is_a_caller_bug() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let mut out = vec![0.0; 5];

        fill_error_spread_divisors(&alleles, &table.view(), &mut out);
    }

    /// **A buffer that is too long is a caller bug too, and only this test says so.** The
    /// check has to be an equality rather than "at least enough": the trailing entries would
    /// never be written, and every genotype past the first would then read a slot the fill
    /// left alone — which is a real number and a wrong one. It is also what makes
    /// [`DivisorTable`]'s bound meaningful, since a longer table admits a genotype index that
    /// should have been out of range.
    #[test]
    #[should_panic(expected = "one entry per (genotype, allele)")]
    fn a_buffer_too_long_is_a_caller_bug_as_well() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let needed = table.view().genotype_count() * 2;
        let mut out = vec![0.0; needed + 1];

        fill_error_spread_divisors(&alleles, &table.view(), &mut out);
    }

    /// A repeat tract's substitution term is a different rate on a different model, so a
    /// divisor table filled for one would be a number with no meaning reaching the wrong row
    /// builder.
    #[test]
    #[should_panic(expected = "the SNP/indel path's")]
    fn a_repeat_tract_is_not_this_paths_locus() {
        use crate::ng::types::Motif;
        let detail = crate::ng::locus_generation::SsrDetail {
            motif: Motif::new(b"AC").expect("AC is a motif"),
            left_flank: b"GGTT"[..].into(),
            right_flank: b"TTGG"[..].into(),
        };
        let alleles = CandidateAlleles::new(b"ACAC"[..].into(), LocusKind::Ssr(detail));
        let table = diploid(1);
        let mut out = vec![0.0; table.view().genotype_count()];

        fill_error_spread_divisors(&alleles, &table.view(), &mut out);
    }

    /// An allele id from another locus, or a genotype index past what the table holds, is
    /// caught when the table is read — which is what [`AlleleId`]'s own documentation promises
    /// and what nothing else here checks.
    #[test]
    #[should_panic(expected = "is past the 2 this locus is called over")]
    fn an_allele_the_locus_does_not_have_is_a_caller_bug() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let out = divisors(&alleles, &table);

        let _ = table_over(&out, &table).at(GenotypeIdx(0), AlleleId(2));
    }

    /// A genotype index from another shape must not quietly resolve: row 3 of a triallelic
    /// diploid table is a genotype this biallelic one does not have, and the two would be
    /// different genotypes even where both tables are long enough to hold the index.
    #[test]
    #[should_panic(expected = "genotype 3 is past the 3 this table was filled for")]
    fn a_genotype_the_table_was_not_filled_for_is_a_caller_bug() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let out = divisors(&alleles, &table);

        // Three genotypes over two alleles: six entries, and no row 3.
        assert_eq!(out.len(), 6);
        let _ = table_over(&out, &table).at(GenotypeIdx(3), AlleleId(0));
    }

    /// The row accessor answers `None` where the pair accessor panics — the same fact, in the
    /// shape a caller that wants to handle it reaches for.
    #[test]
    fn the_row_accessor_refuses_a_genotype_the_table_does_not_hold() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let out = divisors(&alleles, &table);
        let divisors_table = table_over(&out, &table);

        assert!(divisors_table.row_of(GenotypeIdx(2)).is_some());
        assert!(divisors_table.row_of(GenotypeIdx(3)).is_none());
        assert_eq!(divisors_table.allele_count(), 2);
    }

    #[test]
    #[should_panic(expected = "belongs to a different locus")]
    fn a_genotype_table_from_another_locus_is_a_caller_bug() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(3);
        let mut out = vec![0.0; table.view().genotype_count() * 2];

        fill_error_spread_divisors(&alleles, &table.view(), &mut out);
    }

    /// The size the divisor is worth, **taken from a filled table rather than from a literal**:
    /// the difference in what a wrongly-explained read is charged, between an allele that gets
    /// the spread and one that does not, at the same locus and the same genotype.
    ///
    /// It is 1.0986 nats, 4.77 on the Phred scale — the number spec §3.5 uses to argue the
    /// choice matters. Computed from `3.0_f64.ln()` instead, this would be a test of
    /// `f64::ln`: it would pass with the fill deleted and the divisor hardcoded, which is
    /// exactly the shape of test the plan's production differential must also avoid.
    #[test]
    fn the_spread_is_worth_one_point_one_nats_per_wrongly_explained_read() {
        // Allele 1 is one substitution from the reference; allele 2 is a deletion.
        let alleles = locus(&[b"ACGT", b"ACCT", b"AT"]);
        let table = diploid(3);
        let out = divisors(&alleles, &table);
        let hom_ref = genotype_carrying(&table, &[2, 0, 0]);

        let spread = table_over(&out, &table).at(hom_ref, AlleleId(1)).ln()
            - table_over(&out, &table).at(hom_ref, AlleleId(2)).ln();

        assert!(
            (spread - 1.0986).abs() < 5e-5,
            "the spread is {spread} nats"
        );
        assert!(
            (spread * 10.0 / std::f64::consts::LN_10 - 4.77).abs() < 5e-3,
            "the spread is {} Phred",
            spread * 10.0 / std::f64::consts::LN_10
        );
    }
}
