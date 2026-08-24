//! The SNP/indel closed form — how probable one sample's reads are at an ordinary site,
//! given each candidate genotype.
//!
//! Spec §3 in code. This file starts with the one piece the formula needs before it can be
//! written: **where a wrong read's probability goes.**

use crate::ng::calling::{CandidateAlleles, GenotypeIdx, GenotypeTableView};
use crate::ng::locus_generation::LocusKind;
use crate::ng::types::{AlleleId, LogProb};

/// How many bases a misread could have gone to — **three, and it is a physical fact rather
/// than a tuning choice**: three bases to go wrong into, one to come back to.
///
/// Named because it is kept in sync with three other places: spec §3.5, the parameter
/// pre-pass's own noise model, and [`LOG_ERROR_SPREAD`], which is what the row actually reads.
pub const ERROR_SPREAD_BASES: f64 = 3.0;

/// What the row subtracts from a wrong read's charge where the error mass spreads three ways —
/// `ln 3`, about 1.0986 nats or 4.77 Phred.
///
/// **The table stores this and not the 3 it is the logarithm of** (decided at B2). Spec §3.3
/// charges an unexplained observation `q_sum + n·(log scale − log m)`, so `log m` is what the
/// inner loop wants, once per `(observation, genotype)` — and the whole shape of that formula
/// is that **no logarithm is taken inside it**: every logarithm it needs is a property of the
/// allele, the genotype or the read group, and is computed before the loop starts. A table
/// holding `m` would put an `ln` back in, measured at 1.392 ns a term against 0.553 ns.
///
/// It also means *divisor* stops being the right word for what is stored — nobody divides by
/// 1.0986 — which is why the table is [`LogErrorSpreadTable`].
pub const LOG_ERROR_SPREAD: f64 = 1.098_612_288_668_109_7;

/// What it subtracts where the model has nothing to say — **nothing, because the error mass is
/// left unspread.** `ln 1 = 0`, which is the conservative direction and favours the reference.
pub const NO_LOG_ERROR_SPREAD: f64 = 0.0;

/// `log m(a, g)` — how much a wrong read's charge is reduced because the error could have gone
/// several ways (spec §3.5).
///
/// [`LOG_ERROR_SPREAD`] where the observation differs from **every** allele the genotype
/// carries by a substitution at exactly one position, [`NO_LOG_ERROR_SPREAD`] otherwise. In
/// spec §3.5's own terms `m` is `3` and `1`; the table holds their logarithms, because that is
/// what the row's inner loop subtracts and taking the logarithm there would be the one place
/// spec §3.3's closed form has none.
///
/// # Why a spread at all
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
/// counts holds the matching row of spreads at the same offset. [`LogErrorSpreadTable`] is
/// the one way to read it, and it carries the stride so no caller has to supply one.
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
/// fill never touched — a real number and a wrong one. And it is what makes
/// [`LogErrorSpreadTable`]'s own bound meaningful: a longer table admits a genotype index that
/// should have been out of range.
pub fn fill_log_error_spreads(
    alleles: &CandidateAlleles,
    genotypes: &GenotypeTableView<'_>,
    out: &mut [f64],
) {
    // **A repeat tract has no business here.** Its substitution term is a different rate on a
    // different model (spec §4.3), and this spread describes neither — so a table filled for
    // one would be a number with no meaning quietly reaching the wrong row builder.
    assert!(
        matches!(alleles.kind(), LocusKind::Generic),
        "the error spread is the SNP/indel path's; this locus is {:?}",
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
        "the error-spread table needs one entry per (genotype, allele) — {} genotypes × {} alleles \
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
            // spread is never read for an observation the genotype explains anyway.
            let every_carried_is_one_substitution_away = carried_copies
                .iter()
                .enumerate()
                .filter(|&(_, &copies)| copies > 0)
                .all(|(carried_allele, _)| {
                    one_substitution_apart[observed_allele * allele_count + carried_allele]
                });
            out[genotype * allele_count + observed_allele] =
                if every_carried_is_one_substitution_away {
                    LOG_ERROR_SPREAD
                } else {
                    NO_LOG_ERROR_SPREAD
                };
        }
    }
}

/// A filled error-spread table, with the stride it was filled at.
///
/// **The stride travels with the buffer rather than being handed in at each lookup, and that
/// is the whole reason this type exists.** An accessor taking `(values, allele_count, genotype,
/// allele)` cannot check that `allele_count` is the stride the buffer was actually filled at —
/// so reading a three-allele table at a stride of two returns a real spread from the wrong
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
pub struct LogErrorSpreadTable<'a> {
    values: &'a [f64],
    allele_count: usize,
}

impl<'a> LogErrorSpreadTable<'a> {
    /// Wrap a buffer [`fill_log_error_spreads`] filled, against the genotype table it was
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
            "an error-spread table for {} genotypes over {} alleles holds {} entries, not {}",
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

    /// One genotype's whole row of log spreads — one entry per allele, in allele order.
    ///
    /// **The shape the row function wants**, because its inner loop already holds the matching
    /// row of copy counts as a slice and walks the two together.
    ///
    /// `None` where this table holds no such genotype, for the reason [`GenotypeIdx`] gives:
    /// row 4 of a triallelic diploid table and row 4 of a tetraploid one are different
    /// genotypes, so an index from another shape must not quietly resolve.
    #[must_use]
    pub fn log_spreads_for(&self, genotype: GenotypeIdx) -> Option<&'a [f64]> {
        let start = (genotype.get() as usize).checked_mul(self.allele_count)?;
        self.values.get(start..start + self.allele_count)
    }

    /// The log spread for one `(genotype, allele)` pair.
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
        let row = self.log_spreads_for(genotype).unwrap_or_else(|| {
            panic!(
                "genotype {} is past the {} this table was filled for",
                genotype.get(),
                self.values.len() / self.allele_count
            )
        });
        row[allele]
    }

    /// One allele's column of log spreads — its entry for every genotype, in genotype order.
    ///
    /// **The shape the row's inner loop wants**, because for a fixed observation the allele is
    /// fixed and what varies is the genotype. Walking a column checks the allele once and then
    /// strides; calling [`at`](Self::at) per genotype would put a bounds check, a
    /// `checked_mul` and a formatting closure inside a loop whose whole stated shape is a
    /// multiply and an add.
    ///
    /// # Panics
    ///
    /// **In release as well as debug**, on an allele this locus is not called over.
    pub(crate) fn log_spreads_of(&self, allele: AlleleId) -> impl Iterator<Item = f64> + use<'a> {
        let allele = usize::from(allele.get());
        assert!(
            allele < self.allele_count,
            "allele {allele} is past the {} this locus is called over",
            self.allele_count
        );
        self.values[allele..]
            .iter()
            .step_by(self.allele_count)
            .copied()
    }

    /// The same column as [`log_spreads_of`](Self::log_spreads_of), read as `m` rather than as `log m`.
    ///
    /// **This is what the row reads**, and [`log_spreads_of`](Self::log_spreads_of) is what it read
    /// before the contamination mixture arrived. Spec §3.6 evaluates the mixture in
    /// probability space and takes one logarithm of the result (spec §8), so `(1 − c)·ε̄/m`
    /// needs the spread as the number the spec calls `m` — 3 or 1 — and not as its logarithm.
    ///
    /// **So the table's stored form now costs an `exp` per error-side term**, one for every
    /// `(observation, genotype)` pair the genotype cannot explain. The reason it is stored the
    /// other way is recorded and was right when it was taken: spec §3.3's closed form takes no
    /// logarithm inside its loop, so a table of `m` would have put one there (B2, 2026-08-24).
    /// The mixture takes a logarithm inside that loop by specification, which is what makes the
    /// argument no longer hold. Two cheaper shapes are available and neither belongs to this
    /// step: the filler could write `m` (a `log` moves to whoever still wants one, and nothing
    /// in the caller does), or the table could carry both, for one extra buffer of
    /// `genotypes × alleles` — 126 values, 1,008 bytes, at a six-allele diploid locus.
    ///
    /// `exp(0)` is exactly `1.0`, so the spread that does not apply costs nothing in accuracy;
    /// `tests::the_linear_column_is_the_bases_the_log_column_is_the_logarithm_of` pins the
    /// other one against [`ERROR_SPREAD_BASES`].
    ///
    /// # Panics
    ///
    /// **In release as well as debug**, on an allele this locus is not called over.
    pub fn spreads_of(&self, allele: AlleleId) -> impl Iterator<Item = f64> + use<'a> {
        self.log_spreads_of(allele).map(f64::exp)
    }

    /// How many alleles the locus is called over — the table's stride.
    #[must_use]
    pub fn allele_count(&self) -> usize {
        self.allele_count
    }

    /// How many genotypes the table was filled for — the number of rows behind the stride.
    ///
    /// **The row needs this to refuse a table from another ploidy**, and could not ask for it
    /// until now. A row walks a column by striding, and a column shorter than the genotype
    /// count does not fail: it `zip`s short, so the genotypes past the table's end keep
    /// whatever they were seeded with and the row comes back a plausible length. Measured on a
    /// tetraploid handed a diploid's table, the last genotype scored `0.0` against its own
    /// `−0.863` — and it was the winning genotype in the truncated row.
    #[must_use]
    pub fn genotype_count(&self) -> usize {
        self.values.len() / self.allele_count
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

/// One sample's whole `Lg` row at an ordinary site — one log-probability per candidate
/// genotype: spec §3.6's contamination mixture, which is spec §3.3 wherever nothing is
/// contaminated.
///
/// ```text
/// log Lg(g)  =   Σ  n_o · log[ (1 − c_{r(o)}) · own(o | g)  +  c_{r(o)} · q(o) ]
///                o
///
///             +  q_sum_other                                       ← the pooled leftover
/// ```
///
/// **A read either shows something the genotype can produce, in which case it is charged only
/// for which copy it came from, or it does not, in which case it is charged for being wrong** —
/// that is `own(o | g)`, `k_a/P` on the first branch and `ε̄/m` on the second. Beside it sits the
/// chance the read is not this individual's at all: `c` of this read group's reads came from
/// somebody else, and `q(o)` is how often that somebody else shows the allele this observation
/// shows. `q_sum_other` is the same number in every row so it cancels when the caller
/// normalises, and it is kept because the data likelihood also feeds emission and QUAL, where an
/// absolute value is compared between loci.
///
/// # There is no `c == 0` branch, and that is a decision rather than an omission
///
/// With [`ContaminationMixture::none`] the mixture is `1 · own(o | g) + 0`, which is
/// `own(o | g)` bit for bit, and the row computes spec §3.3. It does not compute it *bitwise*:
/// §3.3 charges a wrong read `q_sum + n·(log scale − log m)` in log space where this takes one
/// logarithm of `ε̄/m` in probability space, and `ε̄` is `exp(q_sum/n)` by construction — so the
/// two forms differ by an `exp`/`log` round trip, a few units in the last place.
/// `tests::the_mixture_at_no_contamination_is_the_plain_formula` measures it across the sweep.
/// **Production does keep such a branch and can afford to**: its own two forms differ
/// algebraically, by an extra `(1 − ε)` factor and an allele-count divisor, so it has a real
/// discontinuity to jump rather than a rounding one (spec §3.6).
///
/// # One logarithm per term, and what it costs
///
/// Spec §8 puts exactly one logarithm on each `(observation, genotype)` pair, which is where
/// §3.3's shape — every logarithm hoisted out of the loop — stops applying. Two things are still
/// hoisted, because they can be. The **explained** side depends on the genotype only through its
/// copy count, so its logarithm is taken once per copy count per observation, into a
/// `ploidy + 1` array, rather than once per genotype. The **error** side has one logarithm per
/// term and an `exp` beside it, the second because the spread table stores `log m` where the
/// mixture needs `m` — [`LogErrorSpreadTable::spreads_of`] carries that and what it would take
/// to remove it.
///
/// # Two approximations live in this formula and both have a size
///
/// **A read the genotype explains is not charged for being right.** The exact term is
/// `log(k_a/P) + log(1 − ε)` and the second half is dropped, because recovering it needs the
/// *arithmetic* mean of `ε` over the observation's reads and the evidence carries the geometric
/// one. The omission is at most `n·ε` nats for a genotype explaining `n` reads at rate `ε`, and
/// what matters is the *difference* between two genotypes — 0.002 nats at 20 reads and 1 error
/// in a thousand, 0.75 nats at 300 reads and a poor library at 1 in 20. **Negligible at good
/// chemistry, small but real at bad chemistry and high depth**, and it always favours the
/// genotype that explains more reads.
///
/// **Every read supporting one allele is charged the same, whatever its own quality**, once it
/// is on the error side. That is exactly what `q_sum` encodes and is not an approximation at
/// all: the sum of logs is the log of the product, which is what a labelled-read likelihood
/// wants.
///
/// # What contamination costs the aggregation contract, with its size
///
/// Once `c` is above zero the logarithm sits outside a sum, so the reads pooled into one
/// observation can no longer each carry their own error probability and the geometric mean is
/// substituted for them (spec §1.4). Spec §3.6 measures it: on 20 reads none of which the
/// genotype explains, half at Phred 30 and half at Phred 20, at `c = 0.03` with the contaminant
/// carrying that allele at 1 in 1,000, **0.14 nats — six-tenths of a Phred**. It is zero at
/// `c = 0`, it grows with `c`, and it grows with the contaminant's own frequency for the allele:
/// 1.13 nats at 1 in 100 and 1.89 at 1 in 2. **So it is small where a contaminant allele is rare
/// and not where it is common**, and the aggregation identity spec §12 test 9 pins is a property
/// of the uncontaminated row.
///
/// # What is not here
///
/// **No multinomial coefficient** — production carries one and this does not (spec §3.4), which
/// is a genotype-changing decision the plan measures rather than asserts. **No producer for
/// `q(o)`**: the frequencies arrive as a parameter, and computing them from the batch the sample
/// was sequenced in is the next step (C2). **No partial observations**:
/// [`GenericSampleEvidence::partials`] is not read here, and Milestone D adds the compatibility
/// rule that scores it.
///
/// # Panics
///
/// **In release as well as debug**, on a caller bug (spec §8): an `out` whose length is not the
/// genotype count, a read group with no calibration or no contamination fraction, an observation
/// naming an allele the locus is not called over, or a mixture holding a frequency per allele for
/// a different locus's allele table. Production holds the analogous assertion in release because
/// a scratch array too short for the allele count would otherwise be indexed out of bounds
/// silently.
///
/// [`ContaminationMixture::none`]: super::ContaminationMixture::none
/// [`GenericSampleEvidence::partials`]: super::GenericSampleEvidence::partials
pub fn genotype_log_likelihood_row(
    evidence: &super::GenericSampleEvidence<'_>,
    genotypes: &GenotypeTableView<'_>,
    calibration: &[super::ReadGroupCalibration],
    contamination: super::ContaminationMixture<'_>,
    error_spreads: LogErrorSpreadTable<'_>,
    out: &mut [LogProb],
) {
    let genotype_count = genotypes.genotype_count();
    let allele_count = genotypes.allele_count();
    assert_eq!(
        out.len(),
        genotype_count,
        "a row holds one entry per candidate genotype — {genotype_count}, not {}",
        out.len()
    );
    assert_eq!(
        error_spreads.allele_count(),
        allele_count,
        "the error-spread table is for {} alleles and this genotype table for {allele_count}, so \
         one of them belongs to a different locus",
        error_spreads.allele_count()
    );
    // **The genotype count as well as the stride**, because the inner loop walks a column by
    // striding and `zip` stops at the shorter of the two: a table filled at another ploidy
    // truncates the walk silently, leaving the genotypes past its end holding whatever seeded
    // them. Measured on a tetraploid handed a diploid's table, the last genotype came back
    // `0.0` against its own `−0.863`, and won. (Found by C1's review; the check was missing
    // from B2, which added the stride check beside it.)
    assert_eq!(
        error_spreads.genotype_count(),
        genotype_count,
        "the error-spread table was filled for {} genotypes and this row is over \
         {genotype_count}, so one of them belongs to a different ploidy",
        error_spreads.genotype_count()
    );
    // **Both halves of the mixture are checked against what they index, and both are checked
    // here rather than per observation** — which is what makes `contaminant_frequency_of` and
    // `fraction_of` unable to panic in a run: an allele or a read group an observation names
    // is already known to be one the mixture covers. Lazily, the read-group half would surface
    // only at whichever locus first reached past the end, or never (C1's review).
    assert!(
        contamination.is_absent() || contamination.allele_count() == allele_count,
        "the mixture holds a contaminant frequency for {} alleles and this locus is called over \
         {allele_count}, so one of them belongs to a different locus",
        contamination.allele_count()
    );
    assert!(
        contamination.is_absent() || contamination.read_group_count() == calibration.len(),
        "the mixture holds a fraction for {} read groups and the run supplied {} calibrations, \
         so one of them belongs to a different run",
        contamination.read_group_count(),
        calibration.len()
    );

    // **The pooled leftover seeds every genotype**, which is what makes it cancel: it is the
    // same number in every row. An empty evidence row therefore leaves the row at zero without
    // a branch, because `empty()` carries a zero leftover — the prior decides, which is the
    // right answer rather than a special case (spec §3.3).
    for slot in out.iter_mut() {
        *slot = LogProb(evidence.unmatched_q_sum);
    }

    // `k / P` for every copy count a genotype can carry — the probability a read came from a
    // copy carrying this observation's allele. `k = 0` is never read — that is the error side —
    // and is filled with a value that would be visible if it ever were.
    let copies_of_the_genome = usize::from(genotypes.ploidy().get());
    assert!(
        copies_of_the_genome <= MAX_PLOIDY_COPIES,
        "a sample with {copies_of_the_genome} copies of its genome is past the \
         {MAX_PLOIDY_COPIES} this row builds a copy-share table for"
    );
    let ploidy = f64::from(genotypes.ploidy().get());
    // Filled in its own scope so it is immutable for the rest of the row: it is a property of
    // the ploidy alone, and nothing below may write to it.
    let copy_share = {
        let mut shares = [f64::NAN; MAX_PLOIDY_COPIES + 1];
        for (copies, share) in shares.iter_mut().enumerate().skip(1) {
            *share = copies as f64 / ploidy;
        }
        shares
    };

    // The explained side's logarithms, one per copy count, refilled per observation because both
    // halves of the mixture change with the read group and the allele. Declared here so the row
    // writes into one array rather than standing a fresh one up per observation — the same
    // discipline the caller's own scratch follows, at a scale small enough to live on the stack.
    let mut log_explained_mixture = [f64::NAN; MAX_PLOIDY_COPIES + 1];

    let counts = genotypes.genotype_allele_counts();
    for observation in evidence.supported {
        let allele = usize::from(observation.allele.get());
        assert!(
            allele < allele_count,
            "an observation names allele {allele}, past the {allele_count} this locus is called \
             over"
        );
        let read_group = observation.read_group.get() as usize;
        let scale = calibration.get(read_group).unwrap_or_else(|| {
            panic!(
                "read group {read_group} has no calibration; the run supplied {}",
                calibration.len()
            )
        });
        // Hoisted out of the genotype loop: every one of these is a property of the observation
        // and of the read group it came from, not of the genotype being scored.
        let reads = f64::from(observation.num_reads);
        let contamination_fraction = contamination.fraction_of(observation.read_group);
        let from_this_individual = 1.0 - contamination_fraction;
        // The whole of the mixture's second half for this observation: `c · q(o)` does not
        // depend on the genotype at all, because what the contaminant would have shown is a
        // property of the allele and of the batch, never of what this sample carries.
        // **The read group decides whose frequency this is, not only how much of it.** The
        // batching says which samples ran beside this library, and the contaminant is drawn
        // from those (spec §3.6). Under the default batching every read group reads the same
        // row and this is the cohort frequency.
        let from_somebody_else_carrying_this_allele = contamination_fraction
            * contamination.contaminant_frequency_of(observation.read_group, observation.allele);
        // `(1 − c) · ε̄` — this individual's share of what a read of this observation is charged
        // for being wrong, before the spread divides it. `ε̄` is floored and deliberately not
        // capped; [`ReadGroupCalibration::charged_error`] says why.
        let this_individuals_charged_error =
            from_this_individual * scale.charged_error(observation.q_sum, observation.num_reads);

        // **The explained side varies with the genotype only through its copy count**, so its
        // logarithm is taken here — at most `ploidy` of them per observation — rather than once
        // per genotype. At a six-allele diploid that is 2 logarithms against 21.
        for (log_mixture, share) in log_explained_mixture
            .iter_mut()
            .zip(copy_share)
            .take(copies_of_the_genome + 1)
            .skip(1)
        {
            *log_mixture =
                (from_this_individual * share + from_somebody_else_carrying_this_allele).ln();
        }

        // **Two columns, walked together.** For a fixed observation the allele is fixed, so
        // what the genotype loop needs is this allele's column of copy counts and its column of
        // spreads — both strided by the allele count, both checked once here rather than per
        // term. A bounds-checked accessor per `(observation, genotype)` would put a
        // `checked_mul` and a formatting closure inside the loop.
        let carried_copies = counts[allele..].iter().step_by(allele_count);
        let spreads = error_spreads.spreads_of(observation.allele);

        for ((slot, &copies), spread) in out.iter_mut().zip(carried_copies).zip(spreads) {
            slot.0 += reads
                * if copies > 0 {
                    log_explained_mixture[copies as usize]
                } else {
                    // The one thing that depends on the genotype: how many things this wrong
                    // read could have shown, given what this genotype carries.
                    (this_individuals_charged_error / spread
                        + from_somebody_else_carrying_this_allele)
                        .ln()
                };
        }
    }
}

/// The widest ploidy the row builds a copy-share table for.
///
/// **`Ploidy::try_new` rejects only zero**, so seventeen copies is constructible and would index
/// past the array. The row asserts on it, in release as well as debug, rather than panicking with
/// `index out of bounds` — every other caller bug in this file says in a sentence what went
/// wrong, and this one used to be the exception.
const MAX_PLOIDY_COPIES: usize = 16;

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

    /// Fill the table, so a test reads it through [`LogErrorSpreadTable`] rather than
    /// open-coding the stride.
    fn spreads(alleles: &CandidateAlleles, table: &GenotypeTable) -> Vec<f64> {
        let view = table.view();
        let mut out = vec![0.0; view.genotype_count() * alleles.len()];
        fill_log_error_spreads(alleles, &view, &mut out);
        out
    }

    /// The filled buffer, read through the type that carries its stride.
    fn table_over<'a>(out: &'a [f64], table: &GenotypeTable) -> LogErrorSpreadTable<'a> {
        LogErrorSpreadTable::over(out, &table.view())
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
        let out = spreads(&alleles, &table);
        let hom_ref = genotype_carrying(&table, &[2, 0]);

        assert_eq!(
            table_over(&out, &table).at(hom_ref, AlleleId(1)),
            LOG_ERROR_SPREAD
        );
    }

    /// **The multi-position class.** Two alleles the same length differing at two positions have
    /// no three-way spread: there is no single base that went wrong.
    #[test]
    fn a_read_differing_at_two_positions_gets_no_spread() {
        let alleles = locus(&[b"ACGT", b"ATCT"]);
        let table = diploid(2);
        let out = spreads(&alleles, &table);
        let hom_ref = genotype_carrying(&table, &[2, 0]);

        assert_eq!(
            table_over(&out, &table).at(hom_ref, AlleleId(1)),
            NO_LOG_ERROR_SPREAD
        );
    }

    /// **The indel class.** A deletion has no finite set of things a wrong read could have
    /// shown, so the mass is left unspread — the conservative choice, which favours the
    /// reference.
    #[test]
    fn a_read_carrying_an_indel_gets_no_spread() {
        let alleles = locus(&[b"ACGT", b"AT"]);
        let table = diploid(2);
        let out = spreads(&alleles, &table);
        let hom_ref = genotype_carrying(&table, &[2, 0]);

        assert_eq!(
            table_over(&out, &table).at(hom_ref, AlleleId(1)),
            NO_LOG_ERROR_SPREAD
        );
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
        let out = spreads(&alleles, &table);
        let het_ref_and_deletion = genotype_carrying(&table, &[1, 0, 1]);
        let hom_ref = genotype_carrying(&table, &[2, 0, 0]);

        // Allele 1 is one substitution from allele 0 and a different length from allele 2.
        assert_eq!(
            table_over(&out, &table).at(het_ref_and_deletion, AlleleId(1)),
            NO_LOG_ERROR_SPREAD
        );
        // …and against the reference homozygote alone it does get the spread, so the fixture is
        // not simply one where nothing ever would.
        assert_eq!(
            table_over(&out, &table).at(hom_ref, AlleleId(1)),
            LOG_ERROR_SPREAD
        );
    }

    /// A heterozygote whose two alleles are *both* one substitution from the observation does
    /// get the spread — the other side of the same rule.
    #[test]
    fn a_genotype_whose_every_allele_is_one_substitution_away_gets_the_spread() {
        // Three alleles at one base: the observation `G` is one substitution from both `A` and
        // `C`.
        let alleles = locus(&[b"A", b"C", b"G"]);
        let table = diploid(3);
        let out = spreads(&alleles, &table);
        let het = genotype_carrying(&table, &[1, 1, 0]);

        assert_eq!(
            table_over(&out, &table).at(het, AlleleId(2)),
            LOG_ERROR_SPREAD
        );
    }

    /// An allele the genotype carries never gets the spread, because it differs from itself at
    /// zero positions rather than at one. The spread is not read for such an observation — it
    /// is on the explained side of the formula — but a table that gave it 3.0 would mean the
    /// predicate had been written as *at most* one rather than *exactly* one, which would also
    /// give an identical pair of alleles the spread.
    #[test]
    fn an_allele_the_genotype_carries_is_not_one_substitution_from_itself() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let out = spreads(&alleles, &table);
        let hom_ref = genotype_carrying(&table, &[2, 0]);
        let het = genotype_carrying(&table, &[1, 1]);

        assert_eq!(
            table_over(&out, &table).at(hom_ref, AlleleId(0)),
            NO_LOG_ERROR_SPREAD
        );
        assert_eq!(
            table_over(&out, &table).at(het, AlleleId(0)),
            NO_LOG_ERROR_SPREAD
        );
        assert_eq!(
            table_over(&out, &table).at(het, AlleleId(1)),
            NO_LOG_ERROR_SPREAD
        );
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
        let out = spreads(&alleles, &table);

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
                    LOG_ERROR_SPREAD
                } else {
                    NO_LOG_ERROR_SPREAD
                };
                assert_eq!(
                    table_over(&out, &table)
                        .at(GenotypeIdx(genotype as u32), AlleleId(observed as u16)),
                    expected,
                    "genotype {genotype}, allele {observed}"
                );
                // …and the row accessor agrees with the pair accessor, which is what lets the
                // row function walk a genotype's spreads beside its copy counts.
                assert_eq!(
                    table_over(&out, &table)
                        .log_spreads_for(GenotypeIdx(genotype as u32))
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
        fill_log_error_spreads(&alleles, &view, &mut out);

        let three_ref_one_deletion = genotype_carrying(&tetraploid, &[3, 0, 1]);
        let four_ref = genotype_carrying(&tetraploid, &[4, 0, 0]);

        assert_eq!(
            table_over(&out, &tetraploid).at(three_ref_one_deletion, AlleleId(1)),
            NO_LOG_ERROR_SPREAD
        );
        assert_eq!(
            table_over(&out, &tetraploid).at(four_ref, AlleleId(1)),
            LOG_ERROR_SPREAD
        );
    }

    /// **A locus with only its reference**, which is what `CandidateAlleles::new` produces and
    /// the commonest shape in a genome. One genotype at any ploidy, one entry, and the spread
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
            fill_log_error_spreads(&alleles, &view, &mut out);

            assert_eq!(
                table_over(&out, &table).at(GenotypeIdx(0), AlleleId::REFERENCE),
                NO_LOG_ERROR_SPREAD
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

            fill_log_error_spreads(&alleles, &view, &mut out);

            let spread_table = table_over(&out, &table);
            assert_eq!(spread_table.allele_count(), allele_count);
            // A homozygote's own allele is never one substitution from itself; every other
            // allele is one substitution from it, so every other allele gets the spread.
            let hom_ref_counts: Vec<u32> = std::iter::once(2)
                .chain(std::iter::repeat_n(0, allele_count - 1))
                .collect();
            let hom_ref = genotype_carrying(&table, &hom_ref_counts);
            let row = spread_table
                .log_spreads_for(hom_ref)
                .expect("a genotype the table holds");
            assert_eq!(row[0], NO_LOG_ERROR_SPREAD);
            assert!(row[1..].iter().all(|&spread| spread == LOG_ERROR_SPREAD));
        }
    }

    /// **Every cell is written, whatever the buffer held.** The other tests all pass a freshly
    /// zeroed buffer, where a cell the fill skipped reads as the unspread value, so it
    /// would be caught, but by luck rather than on purpose. Poisoning the buffer with `NaN`
    /// makes a skipped cell fail deliberately, and the buffer is caller scratch reused across
    /// loci, so what it held is the previous locus's answer rather than zero.
    #[test]
    fn every_cell_is_written_over_whatever_the_buffer_held() {
        let alleles = locus(&[b"ACGT", b"ACCT", b"AT"]);
        let table = diploid(3);
        let view = table.view();
        let mut out = vec![f64::NAN; view.genotype_count() * 3];

        fill_log_error_spreads(&alleles, &view, &mut out);

        for (cell, value) in out.iter().enumerate() {
            assert!(
                *value == NO_LOG_ERROR_SPREAD || *value == LOG_ERROR_SPREAD,
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

        fill_log_error_spreads(&alleles, &table.view(), &mut out);
    }

    /// **A buffer that is too long is a caller bug too, and only this test says so.** The
    /// check has to be an equality rather than "at least enough": the trailing entries would
    /// never be written, and every genotype past the first would then read a slot the fill
    /// left alone — which is a real number and a wrong one. It is also what makes
    /// [`LogErrorSpreadTable`]'s bound meaningful, since a longer table admits a genotype
    /// should have been out of range.
    #[test]
    #[should_panic(expected = "one entry per (genotype, allele)")]
    fn a_buffer_too_long_is_a_caller_bug_as_well() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let needed = table.view().genotype_count() * 2;
        let mut out = vec![0.0; needed + 1];

        fill_log_error_spreads(&alleles, &table.view(), &mut out);
    }

    /// A repeat tract's substitution term is a different rate on a different model, so a
    /// spread table filled for one would be a number with no meaning reaching the wrong row
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

        fill_log_error_spreads(&alleles, &table.view(), &mut out);
    }

    /// An allele id from another locus, or a genotype index past what the table holds, is
    /// caught when the table is read — which is what [`AlleleId`]'s own documentation promises
    /// and what nothing else here checks.
    #[test]
    #[should_panic(expected = "is past the 2 this locus is called over")]
    fn an_allele_the_locus_does_not_have_is_a_caller_bug() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let out = spreads(&alleles, &table);

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
        let out = spreads(&alleles, &table);

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
        let out = spreads(&alleles, &table);
        let spread_table = table_over(&out, &table);

        assert!(spread_table.log_spreads_for(GenotypeIdx(2)).is_some());
        assert!(spread_table.log_spreads_for(GenotypeIdx(3)).is_none());
        assert_eq!(spread_table.allele_count(), 2);
    }

    #[test]
    #[should_panic(expected = "belongs to a different locus")]
    fn a_genotype_table_from_another_locus_is_a_caller_bug() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(3);
        let mut out = vec![0.0; table.view().genotype_count() * 2];

        fill_log_error_spreads(&alleles, &table.view(), &mut out);
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
        let out = spreads(&alleles, &table);
        let hom_ref = genotype_carrying(&table, &[2, 0, 0]);

        let spread = table_over(&out, &table).at(hom_ref, AlleleId(1))
            - table_over(&out, &table).at(hom_ref, AlleleId(2));

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

    // ---- B2: the closed form ----

    use super::super::{
        ContaminationMixture, ContaminationView, GenericObservation, GenericSampleEvidence,
        ReadGroupCalibration,
    };
    use crate::ng::parameter_estimation::Provenance;
    use crate::ng::parameter_estimation::joint::contamination::ContaminationSource;
    use crate::ng::types::{BatchId, LogProb, ReadGroupId};

    /// An observation of `num_reads` reads on one allele from one read group, carrying the
    /// summed log error the merge would have folded.
    fn observation(allele: u16, read_group: u32, num_reads: u32, q_sum: f64) -> GenericObservation {
        GenericObservation {
            allele: AlleleId(allele),
            read_group: ReadGroupId(read_group),
            num_reads,
            q_sum,
        }
    }

    /// A calibration whose scale is one, so `log scale` is exactly zero and the arithmetic a
    /// test hand-computes has nothing in it but the terms under examination.
    fn uncalibrated() -> Vec<ReadGroupCalibration> {
        vec![ReadGroupCalibration::defaulted(); 4]
    }

    /// A calibration with a scale that is not one.
    ///
    /// **On its own this catches nothing**, and an earlier version of this comment claimed it
    /// caught a dropped `log scale`: every fixture that used it compared the row against itself
    /// under a fold or a permutation, where the scale cancels. What catches that is
    /// `at_a_scale_of_three_the_calibration_cancels_the_spread` and the production differential,
    /// which now runs at a scale of 2.5.
    fn calibrated(scale: f64) -> Vec<ReadGroupCalibration> {
        vec![
            ReadGroupCalibration {
                scale,
                provenance: Provenance::FittedHere,
            };
            4
        ]
    }

    /// The row on an uncontaminated sample — spec §3.3, which is what most of these fixtures
    /// are about.
    fn row(
        evidence: &GenericSampleEvidence<'_>,
        alleles: &CandidateAlleles,
        table: &GenotypeTable,
        calibration: &[ReadGroupCalibration],
    ) -> Vec<f64> {
        contaminated_row(
            evidence,
            alleles,
            table,
            calibration,
            ContaminationMixture::uncontaminated(),
        )
    }

    fn contaminated_row(
        evidence: &GenericSampleEvidence<'_>,
        alleles: &CandidateAlleles,
        table: &GenotypeTable,
        calibration: &[ReadGroupCalibration],
        contamination: ContaminationMixture<'_>,
    ) -> Vec<f64> {
        let view = table.view();
        let spreads = spreads(alleles, table);
        let mut out = vec![LogProb(f64::NAN); view.genotype_count()];
        genotype_log_likelihood_row(
            evidence,
            &view,
            calibration,
            contamination,
            LogErrorSpreadTable::over(&spreads, &view),
            &mut out,
        );
        out.into_iter().map(LogProb::get).collect()
    }

    /// Enough entries that any fixture's read groups fit; a mixture takes the prefix it needs.
    static ALL_IN_ONE_BATCH: [BatchId; 16] = [BatchId::ONLY; 16];

    /// A mixture under **the default batching** — one batch holding every read group, so every
    /// observation reads the same contaminant frequencies. That is what a run which declares no
    /// batching gets, and it is the shape all but the batch-specific fixtures want.
    fn one_batch<'a>(
        fractions: &'a [ContaminationView],
        frequencies: &'a [f64],
    ) -> ContaminationMixture<'a> {
        ContaminationMixture::new(
            fractions,
            &ALL_IN_ONE_BATCH[..fractions.len()],
            frequencies,
            frequencies.len(),
        )
    }

    /// A contamination fraction for every read group the fixtures use, with the counts that
    /// say it was measured rather than shrugged at.
    fn every_read_group_contaminated_at(fraction: f64) -> Vec<ContaminationView> {
        vec![
            ContaminationView {
                fraction,
                markers_with_reads: 5_000,
                reads_on_markers: 40_000,
                source: ContaminationSource::ThisReadGroupsReads,
            };
            4
        ]
    }

    /// **A hand-computed biallelic diploid case**, every term written out.
    ///
    /// Two reads showing the reference with a summed log error of −6, one read showing the
    /// alternative at −7; the two alleles differ at one position, so a wrong read's charge is
    /// reduced by `ln 3`. With the scale at one:
    ///
    /// - **reference homozygote**: the two reference reads are explained and charged
    ///   `2·ln(2/2) = 0`; the alternative read is an error, charged `−7 + 1·(0 − ln 3)`.
    /// - **heterozygote**: every read is explained, at `ln(1/2)` each — `3·(−ln 2)`.
    /// - **alternative homozygote**: the alternative read is explained at `ln(2/2) = 0`; the two
    ///   reference reads are errors, charged `−6 + 2·(0 − ln 3)`.
    #[test]
    fn a_biallelic_diploid_row_is_what_the_formula_says_term_by_term() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let supported = [observation(0, 0, 2, -6.0), observation(1, 0, 1, -7.0)];
        let evidence = GenericSampleEvidence::new(&supported, 0.0, &[]);

        let scored = row(&evidence, &alleles, &table, &uncalibrated());

        let hom_ref = genotype_carrying(&table, &[2, 0]).get() as usize;
        let het = genotype_carrying(&table, &[1, 1]).get() as usize;
        let hom_alt = genotype_carrying(&table, &[0, 2]).get() as usize;
        let half = 0.5_f64.ln();

        assert!((scored[hom_ref] - (-7.0 - LOG_ERROR_SPREAD)).abs() < 1e-12);
        assert!((scored[het] - 3.0 * half).abs() < 1e-12);
        assert!((scored[hom_alt] - (-6.0 - 2.0 * LOG_ERROR_SPREAD)).abs() < 1e-12);

        // The heterozygote wins here, which is the answer two reads against one should give.
        assert!(scored[het] > scored[hom_ref] && scored[het] > scored[hom_alt]);

        // **And again with a scale that is not one**, because `log scale` is exactly zero at a
        // defaulted calibration — so every assertion above is blind to what the row does with
        // it. At a scale of 2.5 each error-side read is charged `ln 2.5` more.
        let calibrated_scored = row(&evidence, &alleles, &table, &calibrated(2.5));
        let log_scale = 2.5_f64.ln();
        assert!((calibrated_scored[hom_ref] - (-7.0 + log_scale - LOG_ERROR_SPREAD)).abs() < 1e-12);
        assert!(
            (calibrated_scored[hom_alt] - (-6.0 + 2.0 * log_scale - 2.0 * LOG_ERROR_SPREAD)).abs()
                < 1e-12
        );
        // The heterozygote explains every read, so no scale reaches it at all.
        assert_eq!(calibrated_scored[het], scored[het]);
    }

    /// **The pooled leftover is added to every genotype alike**, so it cancels when the caller
    /// normalises and survives for the data likelihood, which compares loci. A row that dropped
    /// it would move no genotype call and every QUAL.
    #[test]
    fn the_pooled_leftover_shifts_every_genotype_by_the_same_amount() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let supported = [observation(0, 0, 2, -6.0), observation(1, 0, 1, -7.0)];

        let without = row(
            &GenericSampleEvidence::new(&supported, 0.0, &[]),
            &alleles,
            &table,
            &uncalibrated(),
        );
        let with = row(
            &GenericSampleEvidence::new(&supported, -4.25, &[]),
            &alleles,
            &table,
            &uncalibrated(),
        );

        for (bare, pooled) in without.iter().zip(&with) {
            assert!((pooled - bare - (-4.25)).abs() < 1e-12);
        }
    }

    /// **A sample that showed nothing scores every genotype at zero, and no branch makes it
    /// so** — an empty sum is zero, so the prior decides alone (spec §3.3).
    #[test]
    fn a_sample_with_no_evidence_scores_every_genotype_at_zero() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);

        let scored = row(
            &GenericSampleEvidence::empty(),
            &alleles,
            &table,
            &uncalibrated(),
        );

        assert!(scored.iter().all(|&value| value == 0.0));
    }

    /// **The aggregation contract, and the reason the formula has this shape** (spec §12 test 9).
    ///
    /// The merge folds every read showing one allele from one read group into a single
    /// observation with a count and a summed log error. That fold must not change the answer —
    /// a likelihood that treated an observation's reads as interchangeable when their qualities
    /// differ would be wrong by an amount nobody measured.
    ///
    /// So: five reads at five *different* error probabilities, scored as five observations of
    /// one read each, and again as the one observation the merge would have built. The reads
    /// are deliberately spread from Phred 40 to Phred 9, because a fixture where they all carry
    /// the same quality cannot tell the two forms apart.
    ///
    /// **The per-read form is built by struct literal rather than through `new`**, because the
    /// constructor requires strictly ascending `(allele, read group)` pairs and five reads on
    /// one allele are five rows sharing a pair. That is the merge's invariant, not the
    /// formula's, and this test is about the formula.
    ///
    /// **The tolerance is [`WORST_AGGREGATION_RELATIVE`] and not zero, though this fixture
    /// happens to agree to the last bit.** Its exactness is a property of five reads at these
    /// qualities, not of the formula: `how_far_the_bitwise_aggregation_claim_reaches` sweeps
    /// wider and finds the two forms apart. Asserting equality here would pin a fixture-specific
    /// accident and break on a change that costs nothing.
    #[test]
    fn pooling_an_observations_reads_does_not_change_the_answer() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let calibration = calibrated(2.5);

        let per_read_errors: Vec<f64> = [40u8, 31, 22, 15, 9]
            .iter()
            .map(|&phred| {
                (-f64::from(phred) / 10.0 * std::f64::consts::LN_10)
                    .exp()
                    .ln()
            })
            .collect();

        let one_at_a_time: Vec<GenericObservation> = per_read_errors
            .iter()
            .map(|&q| observation(1, 0, 1, q))
            .collect();
        let per_read = GenericSampleEvidence {
            supported: &one_at_a_time,
            unmatched_q_sum: 0.0,
            partials: &[],
        };

        let folded = [observation(1, 0, 5, per_read_errors.iter().sum())];
        let aggregate = GenericSampleEvidence::new(&folded, 0.0, &[]);

        let from_reads = row(&per_read, &alleles, &table, &calibration);
        let from_fold = row(&aggregate, &alleles, &table, &calibration);

        for (genotype, (reads, fold)) in from_reads.iter().zip(&from_fold).enumerate() {
            // The hom-alt genotype scores exactly zero here — every read is explained at
            // `ln(2/2)` — so the denominator is floored rather than dividing zero by zero.
            let relative =
                (reads - fold).abs() / reads.abs().max(fold.abs()).max(f64::MIN_POSITIVE);
            assert!(
                relative <= WORST_AGGREGATION_RELATIVE,
                "genotype {genotype}: {reads} from the reads, {fold} from the fold, a relative \
                 {relative} apart"
            );
        }
    }

    /// **How far the bitwise claim actually reaches**, swept rather than asserted from one
    /// fixture — because a single set of qualities agreeing to the last bit could be luck.
    ///
    /// Spec §2.3 says the aggregation identity "*is* bitwise: §3.3's formula sums the same terms
    /// in the same order either way". The sweep below is what that sentence is worth: **144
    /// `(quality set, scale, ploidy, genotype)` combinations** — six quality sets from two to
    /// seven reads and Phred 45 down to 3, three read-group scales including one that is not
    /// one, and ploidy 2 and 4.
    ///
    /// **The two forms disagree, and the specification's sentence was therefore too strong**
    /// (measured 2026-08-24; spec §2.3 is corrected to match). They are not the same sequence of
    /// additions — a read at a time accumulates `q_r + log scale − log m`, where the fold
    /// accumulates `Σ q_r + n·log scale − n·log m` once — so they round differently, and the
    /// sibling test above happens to land on a fixture where they do not. **One fixture agreeing
    /// to the last bit is luck, and this sweep is what tells the two apart.**
    ///
    /// **The bound is relative and not a unit-in-the-last-place count, because the count grows
    /// with depth.** A first version of this sweep stopped at seven reads and put the worst at
    /// two ulps; widening it to the depths this caller commits to — up to 300 reads a position —
    /// takes that past a hundred, because the disagreement is repeated summation and grows with
    /// how much is summed. What does *not* grow is the relative size, which is what
    /// [`WORST_AGGREGATION_RELATIVE`] bounds, and a caller comparing genotypes cares about the
    /// relative one.
    ///
    /// **What the contract actually requires is untouched**, and it is worth separating the two
    /// claims. Spec §2.3's requirement is that pooling must not change the answer *because of
    /// the model* — that no term is a non-linear function of a per-read quality, since `q_sum`
    /// recovers only the geometric mean. The formula's shape guarantees that exactly, and it is
    /// what the test above pins. What is *not* exact is the floating-point summation order.
    const WORST_AGGREGATION_RELATIVE: f64 = 2.0e-14;

    /// How far the aggregation identity holds, measured rather than claimed — see
    /// [`WORST_AGGREGATION_ULPS`].
    #[test]
    fn how_far_the_bitwise_aggregation_claim_reaches() {
        let alleles = locus(&[b"A", b"C"]);

        let mut worst_relative = 0.0_f64;
        let mut worst_ulps = 0u64;
        let mut worst_where = String::new();
        let mut disagreeing = 0usize;
        let mut compared = 0usize;

        for reads in [2usize, 3, 5, 8, 13, 20, 50, 120, 300] {
            for profile in QUALITY_PROFILES {
                let qualities: Vec<u8> = (0..reads).map(|at| profile(at, reads)).collect();
                let errors: Vec<f64> = qualities
                    .iter()
                    .map(|&phred| -f64::from(phred) / 10.0 * std::f64::consts::LN_10)
                    .collect();
                let one_at_a_time: Vec<GenericObservation> =
                    errors.iter().map(|&q| observation(1, 0, 1, q)).collect();
                let folded = [observation(1, 0, reads as u32, errors.iter().sum())];

                for scale in [1.0_f64, 2.5, 0.37] {
                    let calibration = calibrated(scale);
                    for copies in [2u8, 4] {
                        let table = GenotypeTable::build(
                            Ploidy::try_new(copies).expect("a fixture ploidy"),
                            2,
                        );
                        let from_reads = row(
                            &GenericSampleEvidence {
                                supported: &one_at_a_time,
                                unmatched_q_sum: 0.0,
                                partials: &[],
                            },
                            &alleles,
                            &table,
                            &calibration,
                        );
                        let from_fold = row(
                            &GenericSampleEvidence::new(&folded, 0.0, &[]),
                            &alleles,
                            &table,
                            &calibration,
                        );

                        for (genotype, (from_read, from_sum)) in
                            from_reads.iter().zip(&from_fold).enumerate()
                        {
                            let ulps = from_read.to_bits().abs_diff(from_sum.to_bits());
                            let relative = (from_read - from_sum).abs()
                                / from_read.abs().max(from_sum.abs()).max(f64::MIN_POSITIVE);
                            if ulps > 0 {
                                disagreeing += 1;
                            }
                            worst_ulps = worst_ulps.max(ulps);
                            if relative > worst_relative {
                                worst_relative = relative;
                                worst_where = format!(
                                    "{reads} reads, scale {scale}, ploidy {copies}, genotype \
                                     {genotype}: {from_read} from the reads, {from_sum} from the \
                                     fold, {ulps} ulps"
                                );
                            }
                            compared += 1;
                        }
                    }
                }
            }
        }

        // 9 read counts × 5 profiles × 3 scales × (3 diploid + 5 tetraploid genotypes) = 1,080.
        assert_eq!(compared, 1_080, "the sweep covers 1,080 comparisons");

        // **The sweep must actually disagree somewhere**, or it is measuring a number against
        // itself and the bound below means nothing. This is also what pins the row's summation
        // *shape*: an accumulation that summed the read counts as integers and multiplied each
        // logarithm in once at the end is bitwise exact here, so it would drive this to zero —
        // and that would be a different formula from the one spec §2.3 now describes.
        assert!(
            disagreeing > 0,
            "no comparison disagreed, so this sweep pins nothing"
        );

        assert!(
            worst_relative <= WORST_AGGREGATION_RELATIVE,
            "the aggregation identity is off by a relative {worst_relative}, past the \
             {WORST_AGGREGATION_RELATIVE} this sweep measured ({worst_ulps} ulps at worst) — \
             worst at {worst_where}"
        );
    }

    /// Four quality shapes, as a function of a read's position and how many reads there are:
    /// all reads at Phred 30, all at Phred 93 (the highest the read preparation admits), a
    /// 93/1 alternation, and a monotone fall from 45 to 3. The all-equal profiles matter
    /// because they are where repeated addition of one value is least forgiving.
    const QUALITY_PROFILES: [fn(usize, usize) -> u8; 5] = [
        |_, _| 30,
        |_, _| 93,
        |at, _| if at % 2 == 0 { 93 } else { 1 },
        |at, reads| {
            let span = (reads - 1).max(1) as f64;
            (45.0 - 42.0 * (at as f64 / span)).round() as u8
        },
        // **Every read at Phred 1**, which is the only profile here whose *fold* is charged
        // more than an error of a half — 0.794 at a scale of one, 1.99 at 2.5. Every other
        // profile's geometric mean sits far below the ceiling even when single reads do not
        // (the 93/1 alternation folds to 2 × 10⁻⁵), so without this one no fixture in either
        // sweep can tell a capped charge from an uncapped one, and the zero-contamination
        // sweep's claim to catch a reintroduced ceiling was false. Found by C1's review, which
        // reintroduced the cap and watched that sweep pass.
        |_, _| 1,
    ];

    /// **Order independence** (spec §12 test 8): permuting the observations must not move a
    /// genotype's log-likelihood by more than the summation order costs.
    ///
    /// **Not "by a single bit", which is what §12 test 8 says and what this test used to be
    /// named for.** Permuting the observations *is* changing the summation order, so the two
    /// runs round differently — measured on this fixture, one unit in the last place at genotype
    /// 1. The specification is corrected.
    ///
    /// **The property that matters is that the row imposes no order of its own**, which is what
    /// makes a run reproducible at any worker count: the caller always hands it the merge's
    /// order, and the row must not sort, bucket or re-group behind that.
    #[test]
    fn permuting_the_observations_moves_a_row_only_by_what_the_order_costs() {
        let alleles = locus(&[b"A", b"C", b"AT"]);
        let table = diploid(3);
        let calibration = calibrated(1.7);

        let forward = [
            observation(0, 0, 3, -0.3),
            observation(1, 0, 7, -19.25),
            observation(2, 1, 2, -1e-9),
        ];
        let mut reversed = forward;
        reversed.reverse();

        let scored = row(
            &GenericSampleEvidence::new(&forward, -2.5, &[]),
            &alleles,
            &table,
            &calibration,
        );
        let permuted = row(
            &GenericSampleEvidence {
                supported: &reversed,
                unmatched_q_sum: -2.5,
                partials: &[],
            },
            &alleles,
            &table,
            &calibration,
        );

        // The bound is the same relative one the aggregation sweep measured, so the two
        // properties are quoted on one scale rather than two.
        for (genotype, (a, b)) in scored.iter().zip(&permuted).enumerate() {
            let relative = (a - b).abs() / a.abs().max(b.abs()).max(f64::MIN_POSITIVE);
            assert!(
                relative <= WORST_AGGREGATION_RELATIVE,
                "genotype {genotype}: {a} forward, {b} reversed, a relative {relative} apart"
            );
        }
    }

    /// **The production differential** — the whole of Milestone B's claim, and the reason this
    /// step is not merely "a formula that looks right".
    ///
    /// ng's row and production's `standard_log_likelihood` are the same closed form with two
    /// recorded changes: production carries a multinomial coefficient that ng drops (spec §3.4),
    /// and ng divides a wrong read's error mass by three where production divides by nothing
    /// (spec §3.5). Add the coefficient back and take the spread out, and the two must agree —
    /// **every difference attributed, none unexplained.**
    ///
    /// **Three differences, not two**, and the third is why this runs at a scale that is not
    /// one: production has no calibration at all, so ng's `n·log scale` on the error side has to
    /// come back out too. At a defaulted calibration that term is exactly zero and the whole
    /// reconciliation is blind to it — which is how a mutation deleting `log scale` outright
    /// survived every test in this file, at a cost of 372.84 nats on a four-allele fixture.
    ///
    /// **The spread comes out of the table rather than out of a literal `n · ln 3`**, and the
    /// coefficient's `ln(n!)` is written here rather than borrowed from production's table.
    /// Either shortcut would let the comparison agree with itself: the first passes with the
    /// whole error-spread step deleted and the value hardcoded, which is the shape of test B1
    /// shipped and had to repair; the second cancels a wrong table entry on both sides.
    #[test]
    fn ng_and_production_agree_once_the_two_recorded_changes_are_undone() {
        use crate::pileup_record::AlleleSupportStats;
        use crate::var_calling::per_group_merger::standard_log_likelihood;

        /// `ln(n!)`, written here rather than borrowed from production's table.
        ///
        /// Production's own `ln_factorial` would cancel on both sides of the comparison, so a
        /// wrong entry in its table would be invisible to this test — which is the whole thing
        /// a differential is for.
        fn ln_factorial(n: u64) -> f64 {
            (2..=n).map(|i| (i as f64).ln()).sum()
        }

        let alleles = locus(&[b"A", b"C", b"G"]);
        let table = diploid(3);
        let view = table.view();

        // One read group, so ng's per-(allele, read group) rows line up with production's
        // per-allele ones. Awkward counts and qualities on purpose.
        let supported = [
            observation(0, 0, 11, -3.5),
            observation(1, 0, 4, -19.0),
            observation(2, 0, 2, -8.25),
        ];
        let unmatched = -1.75;
        let evidence = GenericSampleEvidence::new(&supported, unmatched, &[]);

        // **A read-group scale that is not one**, because production has none and a row that
        // dropped `log scale` altogether would otherwise reconcile perfectly: at a defaulted
        // calibration `log scale` is exactly zero, so three of this file's hand-computed tests
        // used to be blind to it. Measured worst case for that mutation, on a four-allele
        // fixture at scale 0.37: 372.84 nats, about 1,620 Phred, with every test green.
        let scale = 2.5_f64;
        let log_scale = scale.ln();
        let spread_values = spreads(&alleles, &table);
        let spread_table = LogErrorSpreadTable::over(&spread_values, &view);
        let mut ours = vec![LogProb(f64::NAN); view.genotype_count()];
        genotype_log_likelihood_row(
            &evidence,
            &view,
            &calibrated(scale),
            ContaminationMixture::uncontaminated(),
            spread_table,
            &mut ours,
        );

        let stats: Vec<AlleleSupportStats> = supported
            .iter()
            .map(|o| AlleleSupportStats {
                num_obs: o.num_reads,
                q_sum: o.q_sum,
                fwd: 0,
                placed_left: 0,
                placed_start: 0,
                mapq_sum: 0,
                mapq_sum_sq: 0,
            })
            .collect();
        let other = AlleleSupportStats {
            num_obs: 0,
            q_sum: unmatched,
            fwd: 0,
            placed_left: 0,
            placed_start: 0,
            mapq_sum: 0,
            mapq_sum_sq: 0,
        };

        let mut compared = 0usize;
        for (genotype, ours) in ours.iter().enumerate() {
            let counts = &view.genotype_allele_counts()[genotype * 3..(genotype + 1) * 3];
            // Production takes the genotype as a multiset of allele indices, one per copy.
            let as_copies: Vec<u8> = counts
                .iter()
                .enumerate()
                .flat_map(|(allele, &copies)| std::iter::repeat_n(allele as u8, copies as usize))
                .collect();
            let theirs = standard_log_likelihood(&stats, &other, &as_copies, 3, 2);

            // The coefficient production carries and ng drops, in closed form.
            let carried_reads: u64 = counts
                .iter()
                .enumerate()
                .filter(|&(_, &copies)| copies > 0)
                .map(|(allele, _)| u64::from(supported[allele].num_reads))
                .sum();
            let coefficient = ln_factorial(carried_reads)
                - counts
                    .iter()
                    .enumerate()
                    .filter(|&(_, &copies)| copies > 0)
                    .map(|(allele, _)| ln_factorial(u64::from(supported[allele].num_reads)))
                    .sum::<f64>();

            // The spread ng subtracts and production does not — taken from the table.
            let spread_effect: f64 = supported
                .iter()
                .filter(|o| counts[usize::from(o.allele.get())] == 0)
                .map(|o| {
                    f64::from(o.num_reads) * spread_table.at(GenotypeIdx(genotype as u32), o.allele)
                })
                .sum();

            // Production charges a wrong read its bare `q_sum`; ng charges it
            // `q_sum + n·log scale`. Taking the scale back out is the third and last difference.
            let scale_effect: f64 = supported
                .iter()
                .filter(|o| counts[usize::from(o.allele.get())] == 0)
                .map(|o| f64::from(o.num_reads) * log_scale)
                .sum();

            let reconciled = ours.get() + coefficient + spread_effect - scale_effect;
            assert!(
                (reconciled - theirs).abs() < 1e-9,
                "genotype {genotype} ({as_copies:?}): ng {} + coefficient {coefficient} + \
                 spread {spread_effect} − scale {scale_effect} = {reconciled}, production \
                 {theirs}",
                ours.get()
            );
            compared += 1;
        }
        assert_eq!(
            compared, 6,
            "a diploid over three alleles has six genotypes"
        );
    }

    /// The spread is not zero in that fixture, so the reconciliation is moving a real number.
    ///
    /// **It does not follow that the differential would fail if the fill returned zero
    /// everywhere**, and an earlier version of this comment claimed it would: the differential
    /// adds back the same table value it subtracts, so a uniformly zero table cancels on both
    /// sides. What catches a constant fill is B1's own class tests, and
    /// `a_wrong_read_that_could_not_be_a_misread_takes_no_spread` above.
    #[test]
    fn the_spread_is_not_zero_in_the_production_differentials_fixture() {
        let alleles = locus(&[b"A", b"C", b"G"]);
        let table = diploid(3);
        let out = spreads(&alleles, &table);
        let hom_ref = genotype_carrying(&table, &[2, 0, 0]);

        assert_eq!(
            table_over(&out, &table).at(hom_ref, AlleleId(1)),
            LOG_ERROR_SPREAD
        );
    }

    /// **The named constant is the logarithm of the named count**, so a transcribed digit in one
    /// cannot drift from the other. Nothing else ties them together — [`ERROR_SPREAD_BASES`] is
    /// read by no code at all — and the size test below admits any value within 5 × 10⁻⁵ of
    /// 1.0986, a band a typo in the fifth decimal fits inside.
    #[test]
    fn the_log_spread_is_the_logarithm_of_the_bases_it_is_named_for() {
        assert_eq!(LOG_ERROR_SPREAD, ERROR_SPREAD_BASES.ln());
        assert_eq!(NO_LOG_ERROR_SPREAD, 1.0_f64.ln());
    }

    /// **At a scale of three the calibration and the spread cancel exactly**, so a wrongly
    /// explained read is charged its `q_sum` and nothing else. That is a term-by-term statement
    /// about `log scale` which does not restate `f64::ln`.
    ///
    /// A row that dropped `log scale` charges `q_sum − n·ln 3` here instead. Every fixture with a
    /// scale other than one that came before this compared the row against *itself* under a fold
    /// or a permutation, where the scale cancels — so a review mutation deleting the term outright
    /// survived the whole file.
    #[test]
    fn at_a_scale_of_three_the_calibration_cancels_the_spread() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let supported = [observation(0, 0, 2, -6.0), observation(1, 0, 1, -7.0)];
        let evidence = GenericSampleEvidence::new(&supported, 0.0, &[]);

        let scored = row(&evidence, &alleles, &table, &calibrated(3.0));

        let hom_ref = genotype_carrying(&table, &[2, 0]).get() as usize;
        let hom_alt = genotype_carrying(&table, &[0, 2]).get() as usize;

        assert!(
            (scored[hom_ref] - -7.0).abs() < 1e-12,
            "the wrongly explained read is charged {}, not its own q_sum",
            scored[hom_ref]
        );
        assert!(
            (scored[hom_alt] - -6.0).abs() < 1e-12,
            "the two wrongly explained reads are charged {}, not their own q_sum",
            scored[hom_alt]
        );
    }

    /// **Two read groups a hundredfold apart, and which one the wrong reads came from changes the
    /// genotype called.** Five reference reads from a group at scale 0.01 and six alternative
    /// reads from a group at one: the reference homozygote wins, and a row ignoring `log scale`
    /// calls the alternative homozygote instead.
    #[test]
    fn which_read_group_the_wrong_reads_came_from_changes_the_call() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let supported = [observation(0, 0, 5, -1.0), observation(1, 1, 6, -1.0)];
        let evidence = GenericSampleEvidence::new(&supported, 0.0, &[]);
        let calibration = vec![
            ReadGroupCalibration {
                scale: 0.01,
                provenance: Provenance::FittedHere,
            },
            ReadGroupCalibration::defaulted(),
        ];

        let scored = row(&evidence, &alleles, &table, &calibration);

        let hom_ref = genotype_carrying(&table, &[2, 0]).get() as usize;
        let het = genotype_carrying(&table, &[1, 1]).get() as usize;
        let hom_alt = genotype_carrying(&table, &[0, 2]).get() as usize;

        assert!((scored[hom_ref] - (-1.0 - 6.0 * LOG_ERROR_SPREAD)).abs() < 1e-12);
        assert!((scored[het] - 11.0 * 0.5_f64.ln()).abs() < 1e-12);
        assert!(
            (scored[hom_alt] - (-1.0 + 5.0 * 0.01_f64.ln() - 5.0 * LOG_ERROR_SPREAD)).abs() < 1e-12
        );
        assert!(
            scored[hom_ref] > scored[het] && scored[hom_ref] > scored[hom_alt],
            "the reference homozygote should win: {scored:?}"
        );
    }

    /// **The copy share is the genotype's own ploidy**, term by term at four copies — the one
    /// thing a diploid fixture cannot say, because `ln(k/P)` at `k = P` is zero at every ploidy
    /// and the heterozygote's `ln(1/2)` is all a diploid row constrains.
    ///
    /// A row reading the copy share off a hardcoded ploidy of two returns **positive**
    /// log-probabilities here and reorders the genotypes. The aggregation sweep does not see it,
    /// because it compares a tetraploid row against another tetraploid row.
    #[test]
    fn the_copy_share_is_the_genotypes_own_ploidy() {
        let alleles = locus(&[b"A", b"C"]);
        let table = GenotypeTable::build(Ploidy::try_new(4).expect("four is a ploidy"), 2);
        let supported = [observation(0, 0, 3, -6.0), observation(1, 0, 1, -7.0)];
        let evidence = GenericSampleEvidence::new(&supported, 0.0, &[]);

        let scored = row(&evidence, &alleles, &table, &uncalibrated());
        let at = |copies: &[u32]| scored[genotype_carrying(&table, copies).get() as usize];

        assert!((at(&[4, 0]) - (-7.0 - LOG_ERROR_SPREAD)).abs() < 1e-12);
        assert!((at(&[3, 1]) - (3.0 * 0.75_f64.ln() + 0.25_f64.ln())).abs() < 1e-12);
        assert!((at(&[2, 2]) - 4.0 * 0.5_f64.ln()).abs() < 1e-12);
        assert!((at(&[1, 3]) - (3.0 * 0.25_f64.ln() + 0.75_f64.ln())).abs() < 1e-12);
        assert!((at(&[0, 4]) - (-6.0 - 3.0 * LOG_ERROR_SPREAD)).abs() < 1e-12);
        assert!(
            scored.iter().all(|&value| value < 0.0),
            "a log-probability cannot be positive: {scored:?}"
        );
    }

    /// **A wrongly explained read that could not have been a misread takes no spread, and the row
    /// has to read that out of the table.** Every row fixture before this one used single-base
    /// alleles, where every pair is one substitution apart — so the whole table could be replaced
    /// by the constant `LOG_ERROR_SPREAD` and every test still passed. **That is the gate this
    /// step set for the production differential, still open one level down in the row itself.**
    ///
    /// Here the locus carries a deletion. Against the reference homozygote the substitution allele
    /// takes `ln 3` and the deletion takes nothing; against the deletion homozygote *neither*
    /// wrong read is one substitution away, so the charge is the two `q_sum`s alone.
    #[test]
    fn a_wrong_read_that_could_not_be_a_misread_takes_no_spread() {
        let alleles = locus(&[b"ACGT", b"ACCT", b"AT"]);
        let table = diploid(3);
        let supported = [observation(1, 0, 2, -4.0), observation(2, 0, 3, -5.0)];
        let evidence = GenericSampleEvidence::new(&supported, 0.0, &[]);

        let scored = row(&evidence, &alleles, &table, &uncalibrated());
        let at = |copies: &[u32]| scored[genotype_carrying(&table, copies).get() as usize];

        assert!((at(&[2, 0, 0]) - (-4.0 - 2.0 * LOG_ERROR_SPREAD - 5.0)).abs() < 1e-12);
        assert!(
            (at(&[0, 0, 2]) - -4.0).abs() < 1e-12,
            "a genotype whose every wrong read is an indel away is charged {}, not −4",
            at(&[0, 0, 2])
        );
    }

    /// An error-spread table filled at another locus's stride must not reach the row: its lookups
    /// would return real numbers off the wrong rows.
    #[test]
    #[should_panic(expected = "belongs to a different locus")]
    fn an_error_spread_table_from_another_locus_is_a_caller_bug() {
        let wide = locus(&[b"A", b"C", b"G"]);
        let wide_table = diploid(3);
        let wide_view = wide_table.view();
        let spread_values = spreads(&wide, &wide_table);

        let narrow_table = diploid(2);
        let narrow_view = narrow_table.view();
        let mut out = vec![LogProb(f64::NAN); narrow_view.genotype_count()];

        genotype_log_likelihood_row(
            &GenericSampleEvidence::empty(),
            &narrow_view,
            &uncalibrated(),
            ContaminationMixture::uncontaminated(),
            LogErrorSpreadTable::over(&spread_values, &wide_view),
            &mut out,
        );
    }

    /// An observation naming an allele the locus is not called over reaches the copy counts before
    /// it reaches the table, and would index a neighbouring genotype's count — a real number and a
    /// wrong one. The row says so first.
    #[test]
    #[should_panic(expected = "an observation names allele 3, past the 3")]
    fn an_observation_naming_an_allele_the_locus_lacks_is_a_caller_bug() {
        let alleles = locus(&[b"A", b"C", b"G"]);
        let table = diploid(3);
        let supported = [observation(3, 0, 4, -3.0)];

        let _ = row(
            &GenericSampleEvidence {
                supported: &supported,
                unmatched_q_sum: 0.0,
                partials: &[],
            },
            &alleles,
            &table,
            &uncalibrated(),
        );
    }

    /// A buffer longer than the table it is wrapped against admits a genotype index that should
    /// have been out of range — the case [`LogErrorSpreadTable`]'s own documentation names, and
    /// which only the *fill*'s length check was tested for.
    #[test]
    #[should_panic(expected = "holds 6 entries, not 8")]
    fn a_buffer_longer_than_its_table_cannot_be_wrapped() {
        let table = diploid(2);
        let values = vec![0.0; 8];

        let _ = LogErrorSpreadTable::over(&values, &table.view());
    }

    #[test]
    #[should_panic(expected = "a row holds one entry per candidate genotype")]
    fn a_row_of_the_wrong_length_is_a_caller_bug() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let view = table.view();
        let spread_values = spreads(&alleles, &table);
        let mut out = vec![LogProb(0.0); 2];

        genotype_log_likelihood_row(
            &GenericSampleEvidence::empty(),
            &view,
            &uncalibrated(),
            ContaminationMixture::uncontaminated(),
            LogErrorSpreadTable::over(&spread_values, &view),
            &mut out,
        );
    }

    #[test]
    #[should_panic(expected = "has no calibration")]
    fn an_observation_from_an_uncalibrated_read_group_is_a_caller_bug() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let supported = [observation(1, 9, 1, -7.0)];

        let _ = row(
            &GenericSampleEvidence::new(&supported, 0.0, &[]),
            &alleles,
            &table,
            &uncalibrated(),
        );
    }

    // ---- C1: the contamination mixture ----

    /// How far the mixture at no contamination may sit from spec §3.3's own arithmetic,
    /// **relative**, over the sweep below.
    ///
    /// **A relative bound and not a unit-in-the-last-place count**, for the reason spec §2.3
    /// gives about the aggregation identity: the departure is a rounding of a sum that grows
    /// with depth, so a ulp count grows with it and a relative size does not. Measured worst
    /// across the sweep's 4,440 comparisons: **`7.3 × 10⁻¹⁵`**, at three alleles, 300 reads and
    /// a scale of 2.5 — a log-likelihood of −3.8844873955582386 from the mixture against
    /// −3.884487395558267 from §3.3. This bound is 2.7 times that, the same margin its sibling
    /// [`WORST_AGGREGATION_RELATIVE`] carries. In the units a genotype is decided in, the worst
    /// disagreement is `1 × 10⁻¹³` Phred.
    ///
    /// **The worst case is a small log-likelihood, not a large one**, which is what a relative
    /// bound is for: it arises on the all-Phred-1 profile, where the genotype explains nearly
    /// everything and the row lands near −3.9 rather than near −1,600. The absolute gap there
    /// is `2.8 × 10⁻¹⁴` nats.
    const MIXTURE_AT_ZERO_RELATIVE_TOLERANCE: f64 = 2.0e-14;

    /// **Spec §3.3's closed form, written out again in log space and independently of the
    /// row** — the oracle the mixture is checked against at zero contamination.
    ///
    /// This is the formula the row *used* to be, before spec §3.6's mixture replaced it: every
    /// logarithm outside the observation walk, an explained read charged `n·log(k/P)` and a
    /// wrong one `q_sum + n·(log scale − log m)`. It is deliberately not a call into the row
    /// under a flag — a shared implementation would agree with itself whatever either of them
    /// did.
    ///
    /// **It reads its spreads from [`fill_log_error_spreads`]** rather than from a literal
    /// `ln 3`, so deleting that filler is a compile error here rather than a quiet subtraction
    /// (the gate B1's review put on this file).
    fn plain_log_likelihood_row(
        evidence: &GenericSampleEvidence<'_>,
        alleles: &CandidateAlleles,
        table: &GenotypeTable,
        calibration: &[ReadGroupCalibration],
    ) -> Vec<f64> {
        let view = table.view();
        let allele_count = view.allele_count();
        let ploidy = f64::from(view.ploidy().get());
        let spread_values = spreads(alleles, table);
        let spread_table = LogErrorSpreadTable::over(&spread_values, &view);
        let counts = view.genotype_allele_counts();

        (0..view.genotype_count())
            .map(|genotype| {
                let mut total = evidence.unmatched_q_sum;
                for observation in evidence.supported {
                    let allele = usize::from(observation.allele.get());
                    let copies = counts[genotype * allele_count + allele];
                    let reads = f64::from(observation.num_reads);
                    total += if copies > 0 {
                        reads * (f64::from(copies) / ploidy).ln()
                    } else {
                        let scale = calibration[observation.read_group.get() as usize];
                        let log_spread =
                            spread_table.at(GenotypeIdx(genotype as u32), observation.allele);
                        observation.q_sum + reads * (scale.log_scale() - log_spread)
                    };
                }
                total
            })
            .collect()
    }

    /// **Spec §12's eleventh test: the mixture is the plain formula at zero, and there is no
    /// `c == 0` branch making it so.**
    ///
    /// The row runs its one shipped code path with every fraction at zero; the oracle runs spec
    /// §3.3 in log space. They agree to [`MIXTURE_AT_ZERO_RELATIVE_TOLERANCE`], which is what lets
    /// contamination default on: a clean cohort is untouched by the default.
    ///
    /// **Not bitwise, and the sweep is what shows it.** §3.6's identity is exact algebra — `ε̄`
    /// is `exp(q_sum/n)` by construction, so `n·log ε̄` *is* `q_sum` — but §8 evaluates the
    /// mixture in probability space, which puts an `exp`/`log` round trip between the two forms.
    ///
    /// **What this test fails on, beyond a wrong number.** It fails the moment anyone
    /// reintroduces production's extra `(1 − ε)` factor into `own` (which would move every
    /// explained read) or its allele-count divisor (which would move every wrong one), because
    /// the oracle has neither. **And it fails if the mixture is given a ceiling**, which took
    /// a fifth quality profile to make true: every read at Phred 1 folds to a charged error of
    /// 0.794 at a scale of one and 1.99 at 2.5, so a clamp binds on the row and not on the
    /// log-space oracle. Until that profile was added the claim was false — every other
    /// profile's *fold* sits far below the ceiling however poor its single reads are, and C1's
    /// review reintroduced the cap and watched this sweep pass.
    #[test]
    fn the_mixture_at_no_contamination_is_the_plain_formula() {
        let loci = [
            &[&b"ACGT"[..], &b"ACCT"[..]][..],
            &[&b"ACGT"[..], &b"ACCT"[..], &b"AGGT"[..]][..],
            &[&b"ACGT"[..], &b"ACCT"[..], &b"AT"[..], &b"ACGTT"[..]][..],
        ];

        let mut worst_relative = 0.0_f64;
        let mut worst_where = String::new();
        let mut disagreeing = 0usize;
        let mut compared = 0usize;

        for bases in loci {
            let alleles = locus(bases);
            let allele_count = bases.len();
            for reads in [1usize, 4, 37, 300] {
                for profile in QUALITY_PROFILES {
                    let q_sum: f64 = (0..reads)
                        .map(|at| -f64::from(profile(at, reads)) / 10.0 * std::f64::consts::LN_10)
                        .sum();
                    // One observation on the reference and one on the last alternative, so
                    // every genotype has both an explained and a wrong read to charge.
                    let supported = [
                        observation(0, 0, reads as u32, q_sum),
                        observation((allele_count - 1) as u16, 1, reads as u32, q_sum),
                    ];
                    let evidence = GenericSampleEvidence::new(&supported, -1.75, &[]);

                    for scale in [1.0_f64, 2.5, 0.37] {
                        let calibration = calibrated(scale);
                        for copies in [2u8, 4] {
                            let table = GenotypeTable::build(
                                Ploidy::try_new(copies).expect("a fixture ploidy"),
                                allele_count,
                            );
                            let mixed = row(&evidence, &alleles, &table, &calibration);
                            let plain =
                                plain_log_likelihood_row(&evidence, &alleles, &table, &calibration);

                            for (genotype, (from_mixture, from_plain)) in
                                mixed.iter().zip(&plain).enumerate()
                            {
                                let relative = (from_mixture - from_plain).abs()
                                    / from_mixture
                                        .abs()
                                        .max(from_plain.abs())
                                        .max(f64::MIN_POSITIVE);
                                if relative > 0.0 {
                                    disagreeing += 1;
                                }
                                if relative > worst_relative {
                                    worst_relative = relative;
                                    worst_where = format!(
                                        "{allele_count} alleles, {reads} reads, scale {scale}, \
                                         ploidy {copies}, genotype {genotype}: {from_mixture} \
                                         from the mixture, {from_plain} from §3.3"
                                    );
                                }
                                compared += 1;
                            }
                        }
                    }
                }
            }
        }

        // 4 read counts × 5 profiles × 3 scales, over (3 + 5) genotypes at two alleles,
        // (6 + 15) at three and (10 + 35) at four: 60 × 74 = 4,440.
        assert_eq!(compared, 4_440, "the sweep covers 4,440 comparisons");

        // **The sweep must actually disagree somewhere**, or the row and the oracle are the
        // same arithmetic and this bound means nothing. The explained side *is* bitwise — the
        // mixture multiplies the copy share by one and adds zero — so what has to disagree is
        // the error side, and it is the round trip that makes it.
        assert!(
            disagreeing > 0,
            "no comparison disagreed, so this sweep pins nothing: either the row grew a \
             `c == 0` branch or the oracle is calling it"
        );

        assert!(
            worst_relative <= MIXTURE_AT_ZERO_RELATIVE_TOLERANCE,
            "the mixture at zero contamination is off from §3.3 by a relative \
             {worst_relative}, past the {MIXTURE_AT_ZERO_RELATIVE_TOLERANCE} this sweep measured, on \
             {disagreeing} of {compared} comparisons — worst at {worst_where}"
        );
    }

    /// **One contaminated case, every term written out** (spec §3.6, the plan's second C1 test).
    ///
    /// A diploid at a two-allele locus: two reads showing the reference at a summed log error of
    /// −6, one showing the alternative at −7, and the read groups 3% contaminated. The
    /// contaminating population carries the reference at 999 in 1,000 and the alternative at 1 in
    /// 1,000 — a rare contaminant allele, which is where §3.6 says the mixture costs least and
    /// where the reference homozygote has most to gain.
    ///
    /// For the **reference homozygote**, at a scale of one:
    ///
    /// - the two reference reads are explained, `own = 2/2 = 1`, so each is charged
    ///   `log(0.97 · 1 + 0.03 · 0.999)`;
    /// - the alternative read is not, and the two alleles differ at one position so its error is
    ///   spread three ways: `own = e⁻⁷/3`, charged `log(0.97 · e⁻⁷/3 + 0.03 · 0.001)`.
    ///
    /// **What the second term shows is the whole point of the mixture.** `0.97 · e⁻⁷/3` is
    /// `2.9 × 10⁻⁴` and the contaminant's own contribution is `3.0 × 10⁻⁵` — a tenth as much
    /// again — so a read the genotype cannot explain stops being quite as damning as it was.
    #[test]
    fn a_contaminated_biallelic_diploid_row_is_what_the_mixture_says_term_by_term() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let supported = [observation(0, 0, 2, -6.0), observation(1, 0, 1, -7.0)];
        let evidence = GenericSampleEvidence::new(&supported, 0.0, &[]);

        let fractions = every_read_group_contaminated_at(0.03);
        let frequencies = [0.999, 0.001];
        let scored = contaminated_row(
            &evidence,
            &alleles,
            &table,
            &uncalibrated(),
            one_batch(&fractions, &frequencies),
        );

        let hom_ref = genotype_carrying(&table, &[2, 0]).get() as usize;
        let het = genotype_carrying(&table, &[1, 1]).get() as usize;
        let hom_alt = genotype_carrying(&table, &[0, 2]).get() as usize;

        let (own, other) = (0.97_f64, 0.03_f64);
        let (q_ref, q_alt) = (0.999_f64, 0.001_f64);
        let (error_ref, error_alt) = ((-6.0_f64 / 2.0).exp(), (-7.0_f64).exp());

        // The reference homozygote: both reference reads explained at a full copy share, the
        // alternative read wrong with its error spread three ways.
        let expected_hom_ref = 2.0 * (own * 1.0 + other * q_ref).ln()
            + (own * error_alt / ERROR_SPREAD_BASES + other * q_alt).ln();
        // The heterozygote explains every read at half a copy share, so no error term at all.
        let expected_het =
            2.0 * (own * 0.5 + other * q_ref).ln() + (own * 0.5 + other * q_alt).ln();
        // The alternative homozygote: the alternative read explained, both reference reads
        // wrong — and `q_sum` is −6 over two reads, so each is charged the geometric mean e⁻³.
        let expected_hom_alt = 2.0 * (own * error_ref / ERROR_SPREAD_BASES + other * q_ref).ln()
            + (own * 1.0 + other * q_alt).ln();

        assert!((scored[hom_ref] - expected_hom_ref).abs() < 1e-12);
        assert!((scored[het] - expected_het).abs() < 1e-12);
        assert!((scored[hom_alt] - expected_hom_alt).abs() < 1e-12);
    }

    /// **What contamination is for: it stops one stray read making a heterozygote — and how much
    /// it does that depends almost entirely on how common the stray allele is in the
    /// contaminant.**
    ///
    /// The same fixture as the hand-computed case above — two reference reads and one
    /// alternative read at a 3% fraction — scored against a rising contaminant frequency for the
    /// alternative allele. What is measured is the gap the heterozygote holds over the reference
    /// homozygote, which is what decides the call:
    ///
    /// | the contaminant shows the alternative at | the heterozygote leads by |
    /// |---|---|
    /// | never (no mixture) | 6.019 nats |
    /// | 1 in 1,000 | 5.981 |
    /// | 1 in 100 | 5.376 |
    /// | 1 in 2 | 2.131 |
    ///
    /// **So a rare contaminant allele buys almost nothing and a common one buys 3.9 nats, 17
    /// Phred.** That is the same lever spec §3.6 measures for what the mixture costs the
    /// aggregation contract — 0.14 nats at 1 in 1,000 against 1.89 at 1 in 2 — and it points the
    /// same way, which is worth knowing: the fixture where contamination changes a call most is
    /// the fixture where pooling reads costs most.
    ///
    /// **The monotone fall is the assertion and the sizes are recorded beside it**, because the
    /// sizes belong to this fixture — one wrong read at a summed log error of −7, three copies
    /// against a diploid — and not to the model. §3.6 names the direction to watch: an
    /// overestimated fraction suppresses real heterozygotes by attributing their alternative
    /// reads to the contaminant, and the last row is what that looks like.
    #[test]
    fn a_contaminated_sample_is_pulled_away_from_the_heterozygote() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let supported = [observation(0, 0, 2, -6.0), observation(1, 0, 1, -7.0)];
        let evidence = GenericSampleEvidence::new(&supported, 0.0, &[]);

        let hom_ref = genotype_carrying(&table, &[2, 0]).get() as usize;
        let het = genotype_carrying(&table, &[1, 1]).get() as usize;
        let gap = |scored: &[f64]| scored[het] - scored[hom_ref];

        let clean = gap(&row(&evidence, &alleles, &table, &uncalibrated()));

        let fractions = every_read_group_contaminated_at(0.03);
        let contaminated: Vec<f64> = [0.001_f64, 0.01, 0.5]
            .into_iter()
            .map(|alternative| {
                let frequencies = [1.0 - alternative, alternative];
                gap(&contaminated_row(
                    &evidence,
                    &alleles,
                    &table,
                    &uncalibrated(),
                    one_batch(&fractions, &frequencies),
                ))
            })
            .collect();

        assert!(
            contaminated[0] < clean,
            "contamination should move the reference homozygote toward the heterozygote, and \
             the gap went from {clean} to {} nats",
            contaminated[0]
        );
        assert!(
            contaminated.windows(2).all(|pair| pair[1] < pair[0]),
            "a commoner contaminant allele should close the gap further, and the gaps are \
             {contaminated:?} nats"
        );

        let measured = [clean, contaminated[0], contaminated[1], contaminated[2]];
        for (found, recorded) in measured.into_iter().zip([6.019, 5.981, 5.376, 2.131]) {
            assert!(
                (found - recorded).abs() < 1e-3,
                "the gaps are {measured:?} nats, and the table above says \
                 [6.019, 5.981, 5.376, 2.131]"
            );
        }
    }

    /// **What the mixture does at 300 reads a position, which is the end of the range this
    /// caller commits to and where every other contaminated fixture here stops short.**
    ///
    /// Reads showing the alternative at Phred 60, a 3% fraction, the contaminant carrying that
    /// allele at 1 in 100, scored under a reference homozygote — the genotype that has to
    /// explain them all as errors. Per read, a misread costs `e⁻¹³·⁸/3`, 3 in ten million,
    /// against `0.03 × 0.01`, 3 in ten thousand from the contaminant: **the contaminant route
    /// is 900 times the likelier, so each wrong read is 6.80 nats cheaper than it was.**
    ///
    /// **That is a per-read constant, so it grows linearly with depth** — 20.4 nats at 3 reads
    /// and 2,041 at 300 — and the mechanism is worth stating because it is not the mechanism
    /// the shallow fixtures show. Once `c·q(o)` exceeds `(1 − c)·ε̄/m`, the mixture *floors what
    /// a wrong read can cost*, so evidence against a homozygote stops accumulating at the rate
    /// the read qualities would give. At three reads that is a nudge; at three hundred it is
    /// what decides the call.
    ///
    /// **So an overestimated fraction is worse at depth, not merely as bad**, which is the
    /// direction spec §3.6 names as the one to watch.
    #[test]
    fn the_mixture_floors_what_a_wrong_read_costs_and_that_grows_with_depth() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let hom_ref = genotype_carrying(&table, &[2, 0]).get() as usize;
        let per_read_log_error = -6.0 * std::f64::consts::LN_10;

        let fractions = every_read_group_contaminated_at(0.03);
        let frequencies = [0.99, 0.01];

        let gain_at = |reads: u32| {
            let supported = [observation(
                1,
                0,
                reads,
                per_read_log_error * f64::from(reads),
            )];
            let evidence = GenericSampleEvidence::new(&supported, 0.0, &[]);
            let clean = row(&evidence, &alleles, &table, &uncalibrated());
            let contaminated = contaminated_row(
                &evidence,
                &alleles,
                &table,
                &uncalibrated(),
                one_batch(&fractions, &frequencies),
            );
            contaminated[hom_ref] - clean[hom_ref]
        };

        let (shallow, deep) = (gain_at(3), gain_at(300));

        assert!(
            (shallow - 20.41).abs() < 5e-2 && (deep - 2041.0).abs() < 5.0,
            "the reference homozygote gains {shallow} nats at 3 reads and {deep} at 300"
        );
        // **Linear in depth, because the saving is per read** — the property that makes the
        // shallow fixtures an understatement rather than a smaller version of the same thing.
        assert!(
            (deep / shallow - 100.0).abs() < 1e-6,
            "a hundred times the reads should be a hundred times the gain, and it is \
             {}",
            deep / shallow
        );
    }

    /// A contaminated locus at a ploidy above two, since every other contaminated fixture here
    /// is a diploid and the copy share is the one term that varies with ploidy.
    ///
    /// A tetraploid carrying one copy of the alternative explains a read of it at `1/4`; the
    /// mixture adds the contaminant's own `c · q`, so the reads it explains move too — which a
    /// diploid fixture also shows, but not at a copy share this far from a half.
    #[test]
    fn the_copy_share_and_the_contaminant_mix_at_a_tetraploid_too() {
        let alleles = locus(&[b"A", b"C"]);
        let table = GenotypeTable::build(Ploidy::try_new(4).expect("a fixture ploidy"), 2);
        let one_copy = genotype_carrying(&table, &[3, 1]).get() as usize;
        let supported = [observation(1, 0, 2, -14.0)];
        let evidence = GenericSampleEvidence::new(&supported, 0.0, &[]);

        let fractions = every_read_group_contaminated_at(0.05);
        let frequencies = [0.8, 0.2];
        let scored = contaminated_row(
            &evidence,
            &alleles,
            &table,
            &uncalibrated(),
            one_batch(&fractions, &frequencies),
        );

        // Two reads the genotype explains at a quarter share, each charged
        // `ln(0.95 · 0.25 + 0.05 · 0.2)`.
        let expected = 2.0 * (0.95_f64 * 0.25 + 0.05 * 0.2).ln();
        assert!((scored[one_copy] - expected).abs() < 1e-12);
    }

    /// **The batch decides whose frequency a read is scored against, and two libraries at one
    /// locus get different answers.**
    ///
    /// The same three alternative reads, once from a library that ran beside samples where the
    /// alternative is common (1 in 2) and once from one that ran beside samples where it is
    /// rare (1 in 1,000). Under a reference homozygote the first library's reads are far less
    /// surprising, because a neighbour on that run plausibly carries the allele — **measured,
    /// 12.34 nats, 54 Phred, between two libraries of the same cohort at the same locus with
    /// the same reads and the same contamination fraction.** Per read that is `ln(0.0203 /
    /// 0.000332)`: in the first batch the contaminant route contributes `0.04 × 0.5`, which
    /// swamps the misread's `0.96 × e⁻⁷/3`; in the second it contributes `0.04 × 0.001`, which
    /// does not.
    ///
    /// **A row that ignored the batching would give both the same number**, and that is the
    /// whole of what this step added: the frequency was one vector for the locus and is now one
    /// per batch. It is also why the default batching is safe — with one batch both libraries
    /// land on the same row and the answer is the cohort frequency, which is what
    /// `an_explicit_zero_fraction_is_the_same_as_no_mixture` and the sweep rest on.
    #[test]
    fn two_libraries_of_one_cohort_are_scored_against_their_own_neighbours() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let hom_ref = genotype_carrying(&table, &[2, 0]).get() as usize;

        let fractions = every_read_group_contaminated_at(0.04);
        // Read groups 0 and 1 ran in a batch where the alternative is common; 2 and 3 in one
        // where it is rare. Batch-major, two alleles to a row.
        let batching = [BatchId(0), BatchId(0), BatchId(1), BatchId(1)];
        let frequencies = [0.5, 0.5, 0.999, 0.001];
        let mixture = ContaminationMixture::new(&fractions, &batching, &frequencies, 2);

        let scored_from = |read_group: u32| {
            let supported = [observation(1, read_group, 3, -21.0)];
            contaminated_row(
                &GenericSampleEvidence::new(&supported, 0.0, &[]),
                &alleles,
                &table,
                &uncalibrated(),
                mixture,
            )[hom_ref]
        };

        let beside_carriers = scored_from(0);
        let beside_non_carriers = scored_from(2);

        assert!(
            beside_carriers > beside_non_carriers,
            "the library that ran beside carriers should find three alternative reads less \
             surprising under a reference homozygote — {beside_carriers} against \
             {beside_non_carriers}"
        );
        let gap = beside_carriers - beside_non_carriers;
        assert!(
            (gap - 12.340).abs() < 1e-3,
            "the two libraries differ by {gap} nats"
        );
    }

    /// **One batch is the default and it must lose nothing**: every read group lands on the
    /// same row, so the batching cannot change an answer no matter which read group an
    /// observation came from.
    #[test]
    fn under_one_batch_the_read_group_does_not_change_the_frequency() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let fractions = every_read_group_contaminated_at(0.04);
        let frequencies = [0.7, 0.3];
        let mixture = one_batch(&fractions, &frequencies);

        let scored_from = |read_group: u32| {
            let supported = [observation(1, read_group, 3, -21.0)];
            contaminated_row(
                &GenericSampleEvidence::new(&supported, 0.0, &[]),
                &alleles,
                &table,
                &uncalibrated(),
                mixture,
            )
        };

        // Not merely close: the same row, since the same frequency and the same fraction reach
        // the same arithmetic.
        assert_eq!(scored_from(0), scored_from(3));
    }

    /// A fraction of zero is not a special case, but it must also not be a *different* case
    /// from having no mixture at all — the two ways a caller can say *nothing is contaminated*
    /// have to give one answer.
    #[test]
    fn an_explicit_zero_fraction_is_the_same_as_no_mixture() {
        let alleles = locus(&[b"ACGT", b"ACCT", b"AT"]);
        let table = diploid(3);
        let supported = [
            observation(0, 0, 5, -14.0),
            observation(1, 1, 2, -9.5),
            observation(2, 0, 1, -3.25),
        ];
        let evidence = GenericSampleEvidence::new(&supported, -0.5, &[]);

        let fractions = every_read_group_contaminated_at(0.0);
        let frequencies = [0.6, 0.3, 0.1];
        let explicit = contaminated_row(
            &evidence,
            &alleles,
            &table,
            &calibrated(1.8),
            one_batch(&fractions, &frequencies),
        );
        let absent = row(&evidence, &alleles, &table, &calibrated(1.8));

        assert_eq!(explicit, absent);
    }

    /// Two read groups of one sample, one contaminated and one not: the row must charge each
    /// observation the fraction of the group it came from, and not the sample's.
    ///
    /// **This is the whole reason the fraction takes the read-group grain** (spec §3.6): a
    /// second seedling in the tube contaminates every library alike, an index hop contaminates
    /// one, and only the finer grain can say the second. A row that averaged the two, or read
    /// the first group's fraction for every observation, passes every other test in this file.
    ///
    /// Three alternative reads at a summed log error of −21, against a contaminant carrying that
    /// allele at 1 in 10 and a library 5% contaminated. **The reference homozygote finds them
    /// 8.57 nats — 37 Phred — less surprising in the contaminated library.** A read from the
    /// contaminant explains one at `0.05 × 0.1`, 1 in 200; a misread explains it at `e⁻⁷/3`,
    /// 3 in 10,000 — sixteen times less. With both routes open each read is 17.4 times likelier
    /// than by misreading alone, and the three together are worth `3 × ln 17.4`. The
    /// heterozygote explains them either way and moves 0.12 nats the *other* way, since the
    /// contaminant shows the allele at 1 in 10 where a carried copy shows it at 1 in 2.
    #[test]
    fn each_read_group_is_charged_its_own_fraction() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        // The same wrong read, once from read group 0 and once from read group 1.
        let from_group_0 = [observation(1, 0, 3, -21.0)];
        let from_group_1 = [observation(1, 1, 3, -21.0)];

        let mut fractions = every_read_group_contaminated_at(0.0);
        fractions[1].fraction = 0.05;
        let frequencies = [0.9, 0.1];
        let mixture = one_batch(&fractions, &frequencies);

        let clean_group = contaminated_row(
            &GenericSampleEvidence::new(&from_group_0, 0.0, &[]),
            &alleles,
            &table,
            &uncalibrated(),
            mixture,
        );
        let dirty_group = contaminated_row(
            &GenericSampleEvidence::new(&from_group_1, 0.0, &[]),
            &alleles,
            &table,
            &uncalibrated(),
            mixture,
        );

        let hom_ref = genotype_carrying(&table, &[2, 0]).get() as usize;
        let het = genotype_carrying(&table, &[1, 1]).get() as usize;

        let on_the_error_side = dirty_group[hom_ref] - clean_group[hom_ref];
        let on_the_explained_side = dirty_group[het] - clean_group[het];

        assert!(
            on_the_error_side > 0.0 && on_the_explained_side < 0.0,
            "the contaminated library should find three alternative reads less surprising under \
             a reference homozygote and slightly more so under a heterozygote, and it moved \
             {on_the_error_side} and {on_the_explained_side} nats"
        );
        assert!(
            (on_the_error_side - 8.569).abs() < 1e-3
                && (on_the_explained_side + 0.1225).abs() < 1e-4,
            "the two moves are {on_the_error_side} and {on_the_explained_side} nats, against the \
             8.569 and −0.1225 this test records"
        );
    }

    /// **The allele-range check both column readers rest on**, which had no test until C1's
    /// review disabled it and watched the module stay green.
    ///
    /// It is not a bounds check standing in for a panic that would happen anyway: a stride walk
    /// starting past the first row returns a *shorter* column, not an out-of-range one — three
    /// alleles read at a two-allele stride gives two entries where three genotypes want three —
    /// and the row's `zip` then drops a genotype's term in silence. The row's own assertion
    /// makes this unreachable from `genotype_log_likelihood_row`; the table is public, so the
    /// check has to hold on its own.
    #[test]
    #[should_panic(expected = "past the 2 this locus is called over")]
    fn a_column_for_an_allele_the_table_lacks_is_a_caller_bug() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let values = spreads(&alleles, &table);

        let _ = table_over(&values, &table).spreads_of(AlleleId(2)).count();
    }

    /// The linear column is the bases the log column is the logarithm of.
    ///
    /// **`exp(0)` is exactly one**, so the spread that does not apply costs no accuracy at all;
    /// the one that does comes back within a unit in the last place of
    /// [`ERROR_SPREAD_BASES`], which is what the mixture divides by.
    #[test]
    fn the_linear_column_is_the_bases_the_log_column_is_the_logarithm_of() {
        let alleles = locus(&[b"ACGT", b"ACCT", b"AT"]);
        let table = diploid(3);
        let values = spreads(&alleles, &table);
        let read = table_over(&values, &table);

        for allele in [AlleleId(0), AlleleId(1), AlleleId(2)] {
            for (log_spread, spread) in read.log_spreads_of(allele).zip(read.spreads_of(allele)) {
                if log_spread == NO_LOG_ERROR_SPREAD {
                    assert_eq!(spread, 1.0, "exp(0) has to be exactly one");
                } else {
                    assert!(
                        (spread - ERROR_SPREAD_BASES).abs() <= f64::EPSILON * ERROR_SPREAD_BASES,
                        "the spread came back as {spread}, not {ERROR_SPREAD_BASES}"
                    );
                }
            }
        }
    }

    /// A mixture built for a locus with a different number of alleles would read a real
    /// frequency from the wrong allele, or run off the end — the same failure shape the spread
    /// table's own stride check exists for.
    #[test]
    #[should_panic(expected = "belongs to a different locus")]
    fn a_mixture_from_another_locus_is_a_caller_bug() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let fractions = every_read_group_contaminated_at(0.02);
        let frequencies = [0.5, 0.3, 0.2];

        let _ = contaminated_row(
            &GenericSampleEvidence::empty(),
            &alleles,
            &table,
            &uncalibrated(),
            one_batch(&fractions, &frequencies),
        );
    }

    /// **A spread table filled at another ploidy truncates the walk instead of failing.** The
    /// inner loop strides a column and `zip` stops at the shorter of the two, so a tetraploid
    /// row handed a diploid's table leaves its last genotypes holding what seeded them —
    /// measured at `0.0` against their own `−0.863`, which made the truncated genotype the
    /// winner. The check that refuses it is C1's; the hole was B2's, beside the stride check
    /// that was there.
    #[test]
    #[should_panic(expected = "belongs to a different ploidy")]
    fn an_error_spread_table_from_another_ploidy_is_a_caller_bug() {
        let alleles = locus(&[b"A", b"C"]);
        let diploid_table = diploid(2);
        let diploid_spreads = spreads(&alleles, &diploid_table);

        let tetraploid_table =
            GenotypeTable::build(Ploidy::try_new(4).expect("a fixture ploidy"), 2);
        let tetraploid_view = tetraploid_table.view();
        let mut out = vec![LogProb(f64::NAN); tetraploid_view.genotype_count()];

        genotype_log_likelihood_row(
            &GenericSampleEvidence::empty(),
            &tetraploid_view,
            &uncalibrated(),
            ContaminationMixture::uncontaminated(),
            LogErrorSpreadTable::over(&diploid_spreads, &diploid_table.view()),
            &mut out,
        );
    }

    /// **A mixture that does not cover every read group the run calibrated is refused before
    /// the first observation is scored**, not when one of them happens to name a group past
    /// the end.
    ///
    /// The fixture is what makes the difference visible: four fractions against ten
    /// calibrations, and an observation from read group 1 — which the mixture *does* cover. A
    /// lazy check passes here and would surface only at whichever locus first held a read from
    /// group 4 or beyond, or never. (C1's review; this test asserted the lazy message until
    /// the eager check made it unreachable.)
    #[test]
    #[should_panic(expected = "belongs to a different run")]
    fn a_mixture_that_misses_a_calibrated_read_group_is_a_caller_bug() {
        let alleles = locus(&[b"A", b"C"]);
        let table = diploid(2);
        let supported = [observation(1, 1, 1, -7.0)];
        let fractions = every_read_group_contaminated_at(0.02);
        let frequencies = [0.9, 0.1];
        let calibration = [ReadGroupCalibration::defaulted(); 10];

        let _ = contaminated_row(
            &GenericSampleEvidence::new(&supported, 0.0, &[]),
            &alleles,
            &table,
            &calibration,
            one_batch(&fractions, &frequencies),
        );
    }
}
