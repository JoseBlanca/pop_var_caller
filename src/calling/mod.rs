//! Turning a locus's reads into a call: which alleles are in play, how probable each
//! sample's reads are under each candidate genotype, how likely each genotype is before
//! the reads are looked at, and what comes out. Steps 6 to 9 of the ng caller.
//!
//! **One folder rather than four, and a dependency decided it.** As four sibling
//! folders at the top of `src/`, the three calling architecture documents each had
//! to state by hand that the prior and the likelihood must never import from the
//! inference step — a rule the tree made necessary and could not enforce. Under one
//! folder the shared vocabulary sits here, in `calling/mod.rs`, one level above all
//! four sub-modules, so each imports downward and the rule disappears
//! (`doc/devel/ng/arch/module_layout.md`, principle 1b).
//!
//! Steps 10 and 11 are deliberately elsewhere: phasing and locus filtering read what
//! calling produces rather than taking part in producing it.
//!
//! ## What is here
//!
//! The vocabulary that needs nothing from the sub-modules: [`CandidateAlleles`], the
//! alleles a locus is called over; [`ExpectedAlleleCopies`], the fractional allele
//! counts the loop feeds back to itself; and [`LocusInference`] with
//! [`SampleGenotypeCall`], what a locus produces. Beside them, the three the calling
//! loop takes and gives back: [`LocusEvidence`], one locus's reads per sample and the
//! one place the SNP/indel and repeat-tract paths meet; [`FrozenParameters`], everything
//! the parameter pre-pass fitted, gathered so one borrow crosses the seam; and
//! [`CallingScratch`], every buffer a locus's calling fills, allocated once per worker.
//! Of the four sub-modules — the candidate
//! step, the likelihood, the genotype prior and the inference loop — three are here:
//! [`allele_candidates`], step 6, which narrows the merge's allele table to the sequences
//! worth calling over (`doc/devel/ng/impl_plan/candidate_alleles.md`);
//! [`likelihood`], step 7, how probable this sample's reads are given each genotype
//! (`doc/devel/ng/impl_plan/calling_read_likelihoods.md`); and
//! [`genotype_prior`], step 8, how likely each genotype is before any read is looked at
//! (`doc/devel/ng/impl_plan/calling_prior.md`). **The fourth is [`inference`], step 9, the
//! loop that consumes all three** — it holds the seam every way of handling a cohort crosses
//! and the configuration of the three nested loops, two of which ship switched off
//! (`doc/devel/ng/impl_plan/calling_loop.md`).
//!
//! Beside this file rather than inside any one sub-module: [`genotype_table`], which
//! says which genotypes a locus's alleles make and holds the three flat tables the
//! genotype prior and the read likelihood both read them from. It sits here for the
//! same reason the types above do — two of the four sub-modules consume it, so it
//! belongs one level above them — and its three types are re-exported here, so that
//! all of calling's shared vocabulary is named at one depth. Beside it again, and
//! test-only: `genotype_table_parity`, which checks that table value for value against
//! production's `GenotypeShape`, frozen in `testdata/` at promotion step C5. It is a file of its
//! own rather than a block inside [`genotype_table`]'s own tests so that the table's tests stay
//! free of the comparison.
//!
//! The scalars other steps also name — `AlleleId`, `Phred`, `Genotype` — are not here
//! but in [`crate::types`], with the rest of ng's shared vocabulary.

#[cfg(test)]
mod genotype_table_parity;
#[cfg(test)]
mod loop_parity;
#[cfg(test)]
mod quality_parity;

pub mod allele_candidates;
pub mod evidence_shaping;
pub mod genotype_prior;
pub mod genotype_table;
pub mod inference;
pub mod likelihood;
pub mod parameters_file;
pub mod quality;
pub mod run_parameters;

mod candidate_alleles;
mod frozen_parameters;
mod locus_evidence;

pub use candidate_alleles::{CandidateAlleles, ExpectedAlleleCopies};
pub use frozen_parameters::FrozenParameters;
pub use genotype_table::{GenotypeIdx, GenotypeTable, GenotypeTableView};
pub use likelihood::generic::{
    ERROR_SPREAD_BASES, NO_ERROR_SPREAD, ReferenceBesideLocus, WrongReadCharge,
    WrongReadChargeTable, allele_is_compatible_with_partial, fill_wrong_read_charges,
    partial_row_fits_allele,
};
pub use likelihood::{
    ContaminationMixture, ContaminationView, GenericObservation, GenericRowScratch,
    GenericSampleEvidence, MIN_BASE_ERROR, MIN_CONTAMINANT_FREQUENCY, ReadGroupCalibration,
    ReadGroupParameters, SsrRowScratch, SsrSampleEvidence, fill_batch_allele_copies,
    fill_own_batch_contaminant_row, fill_shared_contaminant_allele_frequencies,
};
use locus_evidence::locus_kind_word;
pub use locus_evidence::{GenericLocusSample, LocusEvidence, SsrLocusShape};

use crate::calling::quality::{ArtifactTestCounts, CountPriorMemo, SiteQualityBuffers};
use crate::locus_generation::LocusKind;
use crate::parameter_estimation::Provenance;
use crate::parameter_estimation::joint::stratum_fits::AlleleLengthVariationRung;
#[cfg(test)]
use crate::types::Ploidy;
use crate::types::{GenomeRegion, Genotype, GenotypeRows, InbreedingF, LogProb, Phred};

/// What every scratch slot holds between [`CallingScratch::prepare_for_locus`] and the pass
/// that writes it: **not a number**, so a slot some pass forgot to fill fails the next
/// arithmetic check it is handed to rather than reaching a genotype as a plausible zero.
///
/// Zero is exactly the mistake to avoid: an expected copy count of zero, a log-prior of
/// zero and a concentration of zero are all legal values a pass could have written, so a
/// zero-filled buffer and a correctly-filled one are indistinguishable. It is the same
/// argument `SsrRowScratch::prepare_emissions` makes for its own fill value.
pub const UNWRITTEN_SCRATCH_VALUE: f64 = f64::NAN;

/// **What building one locus's genotype-likelihood table cost** — the instrument
/// `doc/devel/ng/spec/calling_em_loop.md` §13's test 5 reads.
///
/// **The property it exists to pin is that the expensive part does not grow with the pass
/// count.** A version that recomputed it per pass would give identical genotypes, only slower,
/// so nothing but a counter can tell the two apart (§2, §8).
///
/// # The build has two halves, and this has to say which it is counting
///
/// **With contamination on at an ordinary site, a caller may no longer cache a whole row across
/// iterations**: `q(o)`, the contaminating population's frequency for the allele an observation
/// shows, is the locus's own number and moves with the loop
/// (`doc/devel/ng/spec/read_likelihoods.md` §3.6, corrected 2026-08-24). **A repeat tract's third
/// term does not move** — it is the fit's allele length variation for the tract's stratum, frozen before
/// calling starts (§4.5.1) — so a contaminated tract caches its whole row like an uncontaminated
/// one.
///
/// What it may still cache — and what the first three fields count — is the **emission**: the
/// answer to how one copy of one allele produced one observed sequence. It reads no frequency,
/// so it stays computed once per `(sample, observation, candidate)` per locus (§6.1). What runs
/// again is the per-genotype **assembly** of the row from those emissions, which the last two
/// fields count: one multiply and one add inside a logarithm the row was taking anyway.
///
/// **⚠ "Expensive" is the spec's word for the general case, and on the SNP/indel path with no
/// partial reads it is the wrong way round.** Counted on this module's own three-sample
/// contaminated fixture — 3 rows × 3 observations × 3 alleles, diploid, so 6 genotypes and 9
/// assemblies — the assemblies cost about 486 multiply-adds and 54 logarithms against the
/// emission fill's 9 charged-error calls. What the emission genuinely is expensive for is the
/// per-`(partial read, candidate)` byte comparison, and, on the repeat-tract path, an alignment
/// per `(observation, candidate)`. **What the split buys where the emission is cheap is not
/// arithmetic but the count**: `emission_evaluations` cannot grow with the pass count, so the
/// invariant §13's test 5 pins stays a statement about the model rather than about a fixture.
///
/// **An uncontaminated run assembles once too**, because with no frequency in the formula the
/// assembled row is the same value at every pass — so it keeps exactly the behaviour it had
/// before the split.
///
/// **`emission_evaluations` is accumulated per row fill from the shapes that fill was
/// handed, not counted inside the emission.** The repeat-tract path has an emission seam a
/// counting model could wrap; the SNP/indel path has none — its per-read term is inside the
/// fill. So this counts what the driver *asked for*, which is what catches a rebuild,
/// and does not catch a fill that evaluates the wrong number of terms. Saying which
/// is which is the point of the field names.
#[derive(Clone, Copy, Default, PartialEq, Eq, Debug)]
pub(crate) struct EmissionCost {
    /// **The expensive half.** How many times the locus's frequency-free emissions were
    /// computed for the whole cohort — **one**, whatever the pass count and whether or not the
    /// run fitted a contamination fraction.
    pub(crate) emission_builds: u64,
    /// How many times one sample's emissions were filled — one per row per emission build.
    pub(crate) emission_row_fills: u64,
    /// Summed over those fills, `observations × candidates`: the emission evaluations the
    /// fills were asked for. **`Σ_s`, not a three-way product** — a fixture whose samples
    /// hold equal observation counts is the one shape that hides the difference (§13 test 5).
    pub(crate) emission_evaluations: u64,
    /// **The cheap half.** How many times the whole `rows × genotypes` table was assembled
    /// from those emissions.
    ///
    /// **One wherever the assembled row reads no allele frequency**, whatever the pass count —
    /// which is every uncontaminated locus, **and every repeat tract**, whose third mixture term
    /// is frozen before calling starts (spec §4.5.1). At a **contaminated ordinary site** it is
    /// one for the initialisation pass, one at the head of each pass, and one more against the
    /// settled frequencies before the final pass.
    pub(crate) table_assemblies: u64,
    /// How many rows were assembled, counted one at a time — `rows × table_assemblies` on a
    /// whole table, and less where an assembly stopped short.
    pub(crate) row_assemblies: u64,
}

/// Every buffer one locus's calling fills — **allocated once per worker and reused at every
/// locus**, so calling a locus costs no allocation of its own.
///
/// **The reason is measured, not stylistic.** Production lifted exactly these buffers out of
/// its own iteration after a profile put the allocator's own self-time at about one cycle in
/// six (`src/var_calling/posterior_engine.rs`; `doc/devel/ng/spec/calling_em_loop.md` §8).
///
/// **Candidate selection's buffers live here too** — the per-allele running totals it
/// accumulates and its survivor list — rather than standing alone as a second per-worker
/// object. The same worker selects a locus's alleles and then calls it, so a second
/// allocation would buy nothing (`doc/devel/ng/arch/candidate_alleles.md` §2.4). Nothing
/// about their shape changed on the way in.
///
/// **Sized by [`prepare_for_locus`](Self::prepare_for_locus), one call per locus**, and
/// every buffer comes back holding [`UNWRITTEN_SCRATCH_VALUE`]. Until that call the scratch
/// has no shape, and every accessor refuses rather than handing back an empty slice: a fold
/// over a zero-length buffer writes nothing and sums to a plausible `0.0`, which is the
/// failure the fill value exists to prevent, reached by the door the fill cannot cover. The
/// two row scratches and the candidate-selection buffers own their own sizing.
///
/// **Amortised, not allocation-free.** A buffer grows when a locus is wider than every locus
/// this worker has met and never shrinks back, so a worker ends up holding its widest
/// locus's shape for its lifetime.
///
/// `SsrEmissionScratch` is the repeat-tract emission model's own working memory — the
/// associated `Scratch` of whichever [`SsrEmissionModel`](likelihood::ssr_emission::SsrEmissionModel)
/// the run scores tracts with. **It has no default**, so every construction names the model
/// it is for, the way the `SsrRowScratch` it wraps already does: that seam exists for a
/// bake-off between two emission models, and a defaulted parameter would let a run switch
/// models with no call site changing.
#[derive(Default, Debug)]
pub struct CallingScratch<SsrEmissionScratch> {
    /// Candidate selection's per-allele running totals and its survivor list.
    candidate_selection: allele_candidates::SelectionScratch,
    /// `samples × genotypes`, sample-major: how probable each sample's reads are under each
    /// candidate genotype — the **genotype** likelihood, written `Lg` in
    /// `doc/devel/ng/spec/read_likelihoods.md` §1, as against `Lr`, which is one read
    /// against one allele.
    ///
    /// **Assembled from the emissions below, and how often depends on whether the assembled row
    /// reads an allele frequency.** With no fraction fitted it does not, so the row is the same
    /// value at every pass and is written once per locus. With one fitted, an **ordinary site's**
    /// row holds `q(o)` — the cohort's own frequency for the allele an observation shows — which
    /// moves with the loop, so it is written again at every pass; a **repeat tract's** third term
    /// is the fit's allele length variation for its stratum, frozen before calling starts, so its row is
    /// written once too (`doc/devel/ng/spec/read_likelihoods.md` §3.6, §4.5.1, §6.1). What is
    /// written once in every case is the emissions.
    genotype_likelihoods: Vec<LogProb>,
    /// One entry per genotype: the sample being scored, its log-prior.
    prior_row: Vec<LogProb>,
    /// One entry per genotype: the sample being scored, its posterior over genotypes.
    posterior_row: Vec<f64>,
    /// One entry per allele: the locus's **seed concentration** — one positive number per
    /// allele, read as chromosomes the prior behaves as though it had already seen
    /// ([`genotype_prior`]'s module doc). Built once per locus.
    seed_concentration: Vec<f64>,
    /// One entry per allele: the sample being scored, its own concentration — the seed plus
    /// what the *other* samples showed here.
    sample_concentration: Vec<f64>,
    /// One entry per allele: the genotype prior's own working space, which
    /// `PriorRow::new` borrows mutably while it reads the concentration. **One per allele,
    /// not one per genotype** — the buffer beside it, `prior_row`, is the other, and at a
    /// diploid biallelic locus they are 2 and 3, so each is a legal length for the other.
    /// `PriorRow::new` calls this same buffer `per_allele_scratch`.
    prior_per_allele_workspace: Vec<f64>,
    /// One entry per allele: the **cohort's** expected allele copies as this pass has them,
    /// summed over every sample.
    cohort_expected_copies: Vec<f64>,
    /// One entry per allele: what the previous pass's cohort copies were — the convergence
    /// comparison.
    previous_cohort_expected_copies: Vec<f64>,
    /// `samples × alleles`, sample-major: each **sample's own** expected allele copies. The
    /// prior's leave-one-out term is the cohort's total minus this.
    per_sample_expected_copies: Vec<f64>,
    /// `samples × (ploidy + 1)`, sample-major: the site quality's first step, each sample's
    /// log-likelihood of carrying exactly `c` non-reference copies
    /// (`doc/devel/ng/spec/calling_quality.md` §5.2).
    copy_count_log_likelihoods: Vec<f64>,
    /// The site quality's count axis, and the buffer its fold alternates with. Both are
    /// `ploidy + samples × ploidy + 1` long — **padded by `ploidy` on the left**, so a
    /// copy-count tap reads `[padding − copies]` without a bounds test in the quadratic inner
    /// loop.
    allele_count_distribution: Vec<f64>,
    /// The other half of that alternation.
    allele_count_distribution_next: Vec<f64>,
    /// `samples × ploidy + 1`: the fold's result back in the log domain, then the
    /// unnormalised log-posterior over the cohort's allele count once the prior is applied.
    log_allele_count_distribution: Vec<f64>,
    /// The site quality's Beta-Binomial count prior, kept between loci because it is the same
    /// vector at every locus of a run ([`CountPriorMemo`]).
    count_prior: CountPriorMemo,
    /// One entry per allele: how many reads that allele drew, pooled over the samples the
    /// locus was called on — the walk the artifact summary picks its **primary
    /// alternative** from (`doc/devel/ng/spec/calling_quality.md` §6.3).
    ///
    /// **The one buffer here with no [`UNWRITTEN_SCRATCH_VALUE`] sentinel**, and it is
    /// counts rather than the sentinel's type that decide it: reads are whole, so this is a
    /// `u64` and a `NaN` is not expressible in it. What replaces the sentinel is that the
    /// only thing that reads this buffer is the same call that zeroes and fills it, in that
    /// order, so a stale entry from the previous locus cannot be read. A zero here is a
    /// real answer — an allele no read reached — which is why nothing downstream may treat
    /// it as *unwritten*.
    pooled_allele_reads: Vec<u64>,
    /// `genotypes × alleles`, genotype-major: how a wrong read is charged under each genotype —
    /// the divisor its error mass is spread by and whether the calibration scale applies — the
    /// SNP/indel row's `fill_wrong_read_charges` table.
    ///
    /// **Filled once per locus rather than once per sample**: both are properties of the
    /// candidate sequences and of the genotype being scored, and of nothing a sample showed.
    wrong_read_charges: Vec<WrongReadCharge>,
    /// One entry per **scratch row**: which sample of the run that row holds.
    ///
    /// **The scratch is sized for the samples the locus is called on, not for the run**, so
    /// this is the map back. A sample the candidate step ruled uncallable takes no part in
    /// the loop at all (`doc/devel/ng/spec/calling_em_loop.md` §5.0) and has no row here —
    /// which is what makes the cohort's expected copies, the convergence denominator and the
    /// site quality's count axis all run over the same cohort without any of them being told
    /// to skip anything (§9).
    run_sample_of_each_row: Vec<usize>,
    /// One entry per scratch row: that sample's inbreeding coefficient, in run order.
    ///
    /// A compacted copy of [`FrozenParameters::inbreeding_coefficient_by_sample`], because
    /// the loop walks one coefficient per row and the run's slice has an entry for every
    /// sample including the ones with no row.
    inbreeding_coefficient_by_row: Vec<InbreedingF>,
    /// What the last `prepare_for_locus` reset and the table build has counted since — spec
    /// §13 test 5's instrument.
    emission_cost: EmissionCost,
    /// **One SNP/indel emission cache per scratch row**, and one per row rather than one
    /// reused, because with contamination on the loop reads them again after the next sample
    /// has been scored: the emissions are filled once per locus and the rows are assembled from
    /// them once per pass (`doc/devel/ng/spec/read_likelihoods.md` §6.1).
    ///
    /// **Grown to the locus's row count and never shrunk**, like every other buffer here —
    /// and, like the repeat-tract row's, each one sizes itself per sample, which is why they
    /// are outside [`buffer_fingerprints`](Self::buffer_fingerprints).
    generic_rows: Vec<GenericRowScratch>,
    /// `batches × alleles`, batch-major: how many copies of each allele the samples of each
    /// sequencing batch are expected to carry, this locus, at the loop's current estimate —
    /// the first half of spec §3.6's `q(o)`, and the half that does not depend on which
    /// sample is being scored. Refilled **once per pass**.
    ///
    /// Empty wherever nothing reads it — a run the fit found no contamination in, and **every
    /// repeat tract of any run**, whose mixture reads no per-batch frequency.
    batch_allele_copies: Vec<f64>,
    /// `batches × alleles`, batch-major: the contaminant allele frequencies **every** sample of
    /// the locus shares — the copies above turned into a distribution, nobody left out. Filled
    /// **once per pass**; each sample is scored against these rows for every batch but its own,
    /// whose row is [`own_contaminant_row`](Self::own_contaminant_row).
    ///
    /// Empty wherever nothing reads it — see [`batch_allele_copies`](Self::batch_allele_copies).
    contaminant_allele_frequencies: Vec<f64>,
    /// `alleles`: **one sample's own batch row** — that batch's copies with the sample's own taken
    /// out, as a distribution. Refilled once per sample per pass, because the subtraction is that
    /// sample's.
    own_contaminant_row: Vec<f64>,
    /// The batch [`own_contaminant_row`](Self::own_contaminant_row) was last filled for, so that
    /// a test can see the whole table the last row was scored against.
    own_contaminant_batch: usize,
    /// `run samples × alleles`, sample-major: the expected copies above, scattered from the
    /// scratch rows back onto the **run's** sample axis, with zero at every sample this locus
    /// has no row for.
    ///
    /// **It exists because the batching is the run's and the rows are the locus's.** A batch
    /// whose samples were all ruled uncallable here has no row at all, and summing over rows
    /// would leave that batch's row of the copy table unwritten — which the copy fill refuses,
    /// by name, as a batching that does not describe the run. Scattering onto the run's axis
    /// gives such a batch a row of zeros instead, which is what the M-step already does with an
    /// uncallable sample: it contributes nothing.
    ///
    /// Empty wherever nothing reads it — see [`batch_allele_copies`](Self::batch_allele_copies).
    expected_copies_by_run_sample: Vec<f64>,
    /// How many sequencing batches the two tables above are sized for. **Zero means the
    /// contaminant tables were not prepared**, which is what an uncontaminated locus leaves
    /// behind and what every contaminant accessor refuses.
    batch_count: usize,
    /// How many of the run's samples [`expected_copies_by_run_sample`](Self::expected_copies_by_run_sample)
    /// is sized for. Zero alongside `batch_count`.
    run_sample_count: usize,
    /// The repeat-tract row's own scratch, including the emission model's.
    ssr_row: SsrRowScratch<SsrEmissionScratch>,
    /// **The tract's fitted scoring parameters, one cell per `(read group, candidate)`** —
    /// gathered once per repeat tract and read by every row of it
    /// ([`inference::repeat_tract_parameters`](super::calling::inference::repeat_tract_parameters)).
    ///
    /// **Held here rather than built per locus** for the reason every buffer in this type is:
    /// its seven vectors are cleared and refilled at each tract, so a worker allocates them once
    /// for a whole run. What cannot be held here is the *contexts* built from it — they borrow
    /// these vectors, so a struct owning both would be self-referential — and those are one
    /// allocation per tract, which is the cost the repeat path pays and the SNP/indel path does
    /// not.
    tract_fits: crate::calling::inference::repeat_tract_parameters::TractScoringFits,
    /// How many per-locus slippage re-fit rounds the last locus adopted — zero at the shipped
    /// configuration, and at every SNP/indel locus, where the round is structurally ignored
    /// (`doc/devel/ng/spec/calling_em_loop.md` §5.1).
    ///
    /// **A record for a probe, not an input**: the loop's driver writes it after its rounds
    /// have stopped and nothing on the scoring path reads it. It lives on the scratch because
    /// [`LocusInference`] does not carry it and a test of the round's stopping rule needs
    /// *something* the loop exposes to assert a round count against.
    slippage_refit_rounds: u32,
    /// How many **rows** the buffers above are sized for. Zero means never prepared.
    ///
    /// **Rows, not the run's samples, and the difference is this type's whole shape**: a
    /// locus is prepared for the samples it is *called on*, which is the run's samples minus
    /// the ones the candidate step ruled uncallable
    /// (`doc/devel/ng/spec/calling_em_loop.md` §5.0). `LocusEvidence::sample_count` and
    /// `FrozenParameters::sample_count` both answer the other question, and the two numbers
    /// appear within a few lines of each other in the final pass.
    row_count: usize,
    /// How many genotypes they are sized for.
    genotype_count: usize,
    /// How many alleles they are sized for.
    allele_count: usize,
}

impl<SsrEmissionScratch> CallingScratch<SsrEmissionScratch> {
    /// Size every buffer for one locus and fill it with [`UNWRITTEN_SCRATCH_VALUE`], so
    /// nothing survives from the last one.
    ///
    /// **The shape comes from the locus's own allele table and genotype table, and from
    /// nowhere else.** Taking the two objects rather than two integers is what stops the
    /// allele count and the genotype count being swapped: at a diploid biallelic locus they
    /// are 2 and 3, and every buffer would still be a legal length.
    ///
    /// # Panics
    ///
    /// If the genotype table does not index the locus's alleles. This is the first of the
    /// three caller bugs `doc/devel/ng/spec/calling_em_loop.md` §8 names as assertions, and
    /// this is the one point where the locus's shape is fixed, so it is where the two are
    /// compared. A table built before a discovery round admitted an allele is exactly how
    /// the two come apart; with a table one allele narrow, every per-allele buffer is sized
    /// for the old count and the locus is called over a set that silently excludes the
    /// allele just admitted.
    ///
    /// If `row_count` is zero. Every locus is called on at least one sample, so a locus
    /// prepared for no rows is a run whose sample order went missing rather than a locus
    /// nobody covered.
    ///
    /// If `rows × genotypes` or `rows × alleles` overflows a `usize`.
    /// `row_count` is **one row per sample this locus is called on** — not the run's sample
    /// count, which at a locus with an uncallable sample is larger.
    pub fn prepare_for_locus(
        &mut self,
        row_count: usize,
        alleles: &CandidateAlleles,
        genotypes: &GenotypeTableView<'_>,
    ) {
        assert_eq!(
            genotypes.allele_count(),
            alleles.len(),
            "the genotype table indexes {} alleles and this locus is called over {}, so the \
             table was built for a different allele set — a discovery round that admitted \
             an allele needs its table rebuilt with it",
            genotypes.allele_count(),
            alleles.len()
        );
        assert!(
            row_count > 0,
            "every locus is called on at least one sample, so a locus prepared for no rows is \
             a run whose sample order went missing"
        );
        let genotype_count = genotypes.genotype_count();
        let allele_count = genotypes.allele_count();
        let table_len = row_count.checked_mul(genotype_count).unwrap_or_else(|| {
            panic!(
                "a locus of {row_count} rows over {genotype_count} genotypes needs a \
                 genotype-likelihood table longer than a usize can index"
            )
        });
        let copies_len = row_count.checked_mul(allele_count).unwrap_or_else(|| {
            panic!(
                "a locus of {row_count} rows over {allele_count} alleles needs a \
                 per-sample copies table longer than a usize can index"
            )
        });

        // The last locus's round count must not survive onto a locus that runs none.
        self.slippage_refit_rounds = 0;

        let unwritten = LogProb(UNWRITTEN_SCRATCH_VALUE);
        resize_and_fill(&mut self.genotype_likelihoods, table_len, unwritten);
        resize_and_fill(&mut self.prior_row, genotype_count, unwritten);
        resize_and_fill(
            &mut self.posterior_row,
            genotype_count,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.seed_concentration,
            allele_count,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.sample_concentration,
            allele_count,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.prior_per_allele_workspace,
            allele_count,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.cohort_expected_copies,
            allele_count,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.previous_cohort_expected_copies,
            allele_count,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.per_sample_expected_copies,
            copies_len,
            UNWRITTEN_SCRATCH_VALUE,
        );

        // The site quality's four buffers. **Sized from the genotype table's ploidy**, which
        // is the one place this locus's ploidy is stated, so they cannot disagree with the
        // table the fold walks (`doc/devel/ng/spec/calling_quality.md` §5.2, §9).
        let ploidy = usize::from(genotypes.ploidy().get());
        let copy_count_table_len = row_count.checked_mul(ploidy + 1).unwrap_or_else(|| {
            panic!(
                "a locus of {row_count} rows at ploidy {ploidy} needs a copy-count log-likelihood \
                 table longer than a usize can index"
            )
        });
        let largest_count = row_count.checked_mul(ploidy).unwrap_or_else(|| {
            panic!(
                "a cohort of {row_count} rows at ploidy {ploidy} carries more \
                 chromosomes than a usize can count"
            )
        });
        resize_and_fill(
            &mut self.copy_count_log_likelihoods,
            copy_count_table_len,
            UNWRITTEN_SCRATCH_VALUE,
        );
        let padded_axis = ploidy + largest_count + 1;
        resize_and_fill(
            &mut self.allele_count_distribution,
            padded_axis,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.allele_count_distribution_next,
            padded_axis,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.log_allele_count_distribution,
            largest_count + 1,
            UNWRITTEN_SCRATCH_VALUE,
        );

        // The artifact summary's per-allele read totals. **Zero rather than a sentinel** —
        // see the field's own comment: its only reader zeroes it in the same call.
        resize_and_fill(&mut self.pooled_allele_reads, allele_count, 0);
        let spread_len = genotype_count.checked_mul(allele_count).unwrap_or_else(|| {
            panic!(
                "a locus of {genotype_count} genotypes over {allele_count} alleles needs an \
                 wrong-read charge table longer than a usize can index"
            )
        });
        // **No sentinel here, unlike the `f64` buffers**: an enum has no `NaN`. A cell the fill
        // skipped reads as `OtherwiseDifferent` — main's conservative charge — and the fill
        // refuses a table it would not cover entirely, so none is skipped.
        resize_and_fill(
            &mut self.wrong_read_charges,
            spread_len,
            WrongReadCharge::OtherwiseDifferent,
        );
        // **The two per-row maps are cleared, not resized.** Their length is what says how
        // many rows the caller has claimed so far, and it claims them one at a time after this
        // call — so a resize here would be a length nothing yet means.
        self.run_sample_of_each_row.clear();
        self.inbreeding_coefficient_by_row.clear();
        self.emission_cost = EmissionCost::default();

        // **One emission cache per row**, grown and never shrunk. Each sizes itself when it is
        // filled, so nothing here says how wide a sample is.
        if self.generic_rows.len() < row_count {
            self.generic_rows
                .resize_with(row_count, GenericRowScratch::default);
        }

        // **The contaminant tables are un-sized here and re-sized by the driver**, and only on a
        // run that fitted a fraction. Clearing rather than leaving them is what stops a locus
        // called without contamination from reading the last contaminated locus's frequencies:
        // every accessor refuses a zero batch count, and `resize` alone would have left a table
        // of the right length holding another locus's numbers.
        self.batch_count = 0;
        self.run_sample_count = 0;
        self.batch_allele_copies.clear();
        self.contaminant_allele_frequencies.clear();
        self.own_contaminant_row.clear();
        self.expected_copies_by_run_sample.clear();

        self.row_count = row_count;
        self.genotype_count = genotype_count;
        self.allele_count = allele_count;
    }

    /// How many **rows** the buffers are currently sized for — one per sample this locus is
    /// called on, which at a locus with an uncallable sample is fewer than the run's samples.
    #[inline]
    #[must_use]
    pub fn row_count(&self) -> usize {
        self.row_count
    }

    #[cfg(test)]
    /// How many genotypes the buffers are currently sized for.
    #[inline]
    #[must_use]
    pub fn genotype_count(&self) -> usize {
        self.genotype_count
    }

    #[cfg(test)]
    /// How many alleles the buffers are currently sized for.
    #[inline]
    #[must_use]
    pub fn allele_count(&self) -> usize {
        self.allele_count
    }

    /// One sample's row of the genotype-likelihood table — one entry per candidate genotype.
    ///
    /// **The one spelling of the index.** The table is flat and sample-major, so a caller
    /// slicing it itself would, on a slip, read a different sample's likelihoods and score
    /// that sample's reads onto this one's genotype with nothing failing.
    ///
    /// Until this locus's pass writes it, every entry is [`UNWRITTEN_SCRATCH_VALUE`].
    ///
    /// # Panics
    ///
    /// On an unprepared scratch, or a sample past the count
    /// [`prepare_for_locus`](Self::prepare_for_locus) was given.
    #[inline]
    #[must_use]
    pub fn sample_genotype_likelihoods(&self, sample: usize) -> &[LogProb] {
        let range = self.genotype_row_range(sample);
        &self.genotype_likelihoods[range]
    }

    /// One sample's row of the genotype-likelihood table, to fill.
    ///
    /// # Panics
    ///
    /// As [`sample_genotype_likelihoods`](Self::sample_genotype_likelihoods).
    #[inline]
    pub fn sample_genotype_likelihoods_mut(&mut self, sample: usize) -> &mut [LogProb] {
        let range = self.genotype_row_range(sample);
        &mut self.genotype_likelihoods[range]
    }

    /// One sample's own expected allele copies — one entry per allele. The cohort's sum is
    /// [`cohort_expected_copies`](Self::cohort_expected_copies), and the prior's
    /// leave-one-out term is the difference of the two.
    ///
    /// Until this locus's pass writes it, every entry is [`UNWRITTEN_SCRATCH_VALUE`].
    ///
    /// # Panics
    ///
    /// On an unprepared scratch, or a sample past the count
    /// [`prepare_for_locus`](Self::prepare_for_locus) was given.
    #[inline]
    #[must_use]
    pub fn sample_expected_copies(&self, sample: usize) -> &[f64] {
        let range = self.allele_row_range(sample);
        &self.per_sample_expected_copies[range]
    }

    /// One sample's own expected allele copies, to fill.
    ///
    /// # Panics
    ///
    /// As [`sample_expected_copies`](Self::sample_expected_copies).
    #[inline]
    pub fn sample_expected_copies_mut(&mut self, sample: usize) -> &mut [f64] {
        let range = self.allele_row_range(sample);
        &mut self.per_sample_expected_copies[range]
    }

    /// The **cohort's** expected allele copies as this pass has them, summed over every
    /// sample. Until this pass writes it, every entry is [`UNWRITTEN_SCRATCH_VALUE`].
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub fn cohort_expected_copies(&self) -> &[f64] {
        self.assert_prepared();
        &self.cohort_expected_copies
    }

    #[cfg(test)]
    /// The cohort's expected allele copies, to fill.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    pub fn cohort_expected_copies_mut(&mut self) -> &mut [f64] {
        self.assert_prepared();
        &mut self.cohort_expected_copies
    }

    /// What the previous pass's cohort expected copies were — the convergence comparison.
    /// Before the first pass has advanced, every entry is [`UNWRITTEN_SCRATCH_VALUE`], and
    /// every comparison against it is therefore false.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub fn previous_cohort_expected_copies(&self) -> &[f64] {
        self.assert_prepared();
        &self.previous_cohort_expected_copies
    }

    /// Make this pass's cohort copies the previous pass's, and hand back a buffer to write
    /// the next pass into.
    ///
    /// **A swap rather than a copy**, so there is one spelling of *which buffer is now
    /// which*. The returned buffer arrives holding [`UNWRITTEN_SCRATCH_VALUE`] in every
    /// entry: a pass that leaves an allele unwritten fails the next check it reaches,
    /// rather than reading the pass-before-last's value and letting the convergence test
    /// compare a number no pass wrote.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    pub fn advance_cohort_expected_copies(&mut self) -> &mut [f64] {
        self.assert_prepared();
        std::mem::swap(
            &mut self.cohort_expected_copies,
            &mut self.previous_cohort_expected_copies,
        );
        self.cohort_expected_copies.fill(UNWRITTEN_SCRATCH_VALUE);
        &mut self.cohort_expected_copies
    }

    /// The locus's seed concentration. Until this locus fills it, every entry is
    /// [`UNWRITTEN_SCRATCH_VALUE`].
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub fn seed_concentration(&self) -> &[f64] {
        self.assert_prepared();
        &self.seed_concentration
    }

    /// The locus's seed concentration, to fill once per locus.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    pub fn seed_concentration_mut(&mut self) -> &mut [f64] {
        self.assert_prepared();
        &mut self.seed_concentration
    }

    #[cfg(test)]
    /// The sample being scored, its own concentration.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub fn sample_concentration(&self) -> &[f64] {
        self.assert_prepared();
        &self.sample_concentration
    }

    #[cfg(test)]
    /// The sample being scored, its own concentration, to fill per sample per pass.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    pub fn sample_concentration_mut(&mut self) -> &mut [f64] {
        self.assert_prepared();
        &mut self.sample_concentration
    }

    #[cfg(test)]
    /// The sample being scored, its log-prior over every candidate genotype.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub fn prior_row(&self) -> &[LogProb] {
        self.assert_prepared();
        &self.prior_row
    }

    #[cfg(test)]
    /// The sample's log-prior row, to fill.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    pub fn prior_row_mut(&mut self) -> &mut [LogProb] {
        self.assert_prepared();
        &mut self.prior_row
    }

    /// The sample being scored, its posterior over every candidate genotype.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub fn posterior_row(&self) -> &[f64] {
        self.assert_prepared();
        &self.posterior_row
    }

    #[cfg(test)]
    /// The sample's posterior row, to fill.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    pub fn posterior_row_mut(&mut self) -> &mut [f64] {
        self.assert_prepared();
        &mut self.posterior_row
    }

    #[cfg(test)]
    /// The genotype prior's per-allele working space.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub fn prior_per_allele_workspace(&self) -> &[f64] {
        self.assert_prepared();
        &self.prior_per_allele_workspace
    }

    #[cfg(test)]
    /// The genotype prior's per-allele working space, to lend it.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    pub fn prior_per_allele_workspace_mut(&mut self) -> &mut [f64] {
        self.assert_prepared();
        &mut self.prior_per_allele_workspace
    }

    /// Candidate selection's buffers. **Not sized by
    /// [`prepare_for_locus`](Self::prepare_for_locus)** — it owns its own, and selection
    /// runs before the locus's shape is known.
    #[inline]
    pub fn candidate_selection_mut(&mut self) -> &mut allele_candidates::SelectionScratch {
        &mut self.candidate_selection
    }

    /// The three buffers the SNP/indel row builder needs, for one scratch row.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch, and on a row past the ones it was prepared for.
    #[must_use]
    pub(crate) fn generic_row_buffers_mut(&mut self, row: usize) -> GenericRowBuffers<'_> {
        let range = self.genotype_row_range(row);
        GenericRowBuffers {
            wrong_read_charges: &self.wrong_read_charges,
            own_contaminant_row: &self.own_contaminant_row,
            row_scratch: &mut self.generic_rows[row],
            genotype_likelihoods: &mut self.genotype_likelihoods[range],
        }
    }

    /// The SNP/indel row's per-allele wrong-read charge table, to fill once per locus.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    pub(crate) fn wrong_read_charges_mut(&mut self) -> &mut [WrongReadCharge] {
        self.assert_prepared();
        &mut self.wrong_read_charges
    }

    /// Claim a scratch row for one run sample, with its inbreeding coefficient.
    ///
    /// Called once per sample the locus is **called on**, in the run's sample order, and
    /// **after** [`prepare_for_locus`](Self::prepare_for_locus), which clears the previous
    /// locus's map. Claiming first throws the claims away and leaves the locus with none,
    /// which the map's first reader then refuses by name — measured, "0 of them were
    /// claimed", at the first locus.
    pub(crate) fn claim_row_for(&mut self, run_sample: usize, inbreeding: InbreedingF) {
        self.run_sample_of_each_row.push(run_sample);
        self.inbreeding_coefficient_by_row.push(inbreeding);
    }

    /// Which run sample each scratch row holds, in row order.
    ///
    /// # Panics
    ///
    /// If the claimed rows are not one per prepared row. **The check is here rather than at
    /// the claim**, because a caller that prepared for one shape and claimed for another has
    /// a map whose every entry is a different sample's — and the symptom would be one
    /// sample's reads scored against another's likelihoods rather than an index out of range.
    #[inline]
    #[must_use]
    pub(crate) fn run_sample_of_each_row(&self) -> &[usize] {
        self.assert_rows_claimed();
        &self.run_sample_of_each_row
    }

    /// The claimed rows' inbreeding coefficients, one per row — what the loop walks.
    ///
    /// # Panics
    ///
    /// As [`run_sample_of_each_row`](Self::run_sample_of_each_row).
    #[inline]
    #[must_use]
    pub(crate) fn inbreeding_coefficient_by_row(&self) -> &[InbreedingF] {
        self.assert_rows_claimed();
        &self.inbreeding_coefficient_by_row
    }

    /// Refuse a scratch whose per-row map is not one entry per prepared row.
    #[inline]
    fn assert_rows_claimed(&self) {
        self.assert_prepared();
        assert_eq!(
            self.run_sample_of_each_row.len(),
            self.row_count,
            "this scratch was prepared for {} rows and {} of them were claimed, so the map \
             from rows back to the run's samples describes a different locus",
            self.row_count,
            self.run_sample_of_each_row.len()
        );
    }

    /// What this locus's table build has cost so far — the instrument spec §13's test 5 reads.
    ///
    /// **Nothing in a run reads it**, and that is what it is for: the cost it records is
    /// invisible in the output, because emissions recomputed every pass give the same
    /// genotypes.
    /// Step D2 of the calling loop's plan is what asserts on it.
    #[inline]
    #[must_use]
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "an instrument: what it measures — a table built once rather than once \
                      per pass — is invisible in a run's output, so its readers are step D2's \
                      tests and the bakeoffs plan's benches"
        )
    )]
    pub(crate) fn emission_cost(&self) -> EmissionCost {
        self.emission_cost
    }

    /// Record one sample's emission fill over `observations` observations and `candidates`
    /// candidates.
    pub(crate) fn charge_emission_row_fill(&mut self, observations: usize, candidates: usize) {
        self.emission_cost.emission_row_fills += 1;
        self.emission_cost.emission_evaluations += (observations * candidates) as u64;
    }

    /// Record that the locus's frequency-free emissions were computed once more.
    pub(crate) fn charge_emission_build(&mut self) {
        self.emission_cost.emission_builds += 1;
    }

    /// Record that a fresh assembly of the whole `rows × genotypes` table has begun.
    ///
    /// **The rows are charged one at a time by [`charge_row_assembly`](Self::charge_row_assembly)
    /// rather than added up here**, so that the two counts can disagree. Charged from the row
    /// count, `row_assemblies` is `table_assemblies × rows` by construction and can only repeat
    /// what the other says — measured: an assembly that stops one row short of the table then
    /// leaves both counts untouched, and only the genotypes it did not write say so.
    pub(crate) fn charge_table_assembly(&mut self) {
        self.emission_cost.table_assemblies += 1;
    }

    /// Record one row of the current assembly.
    pub(crate) fn charge_row_assembly(&mut self) {
        self.emission_cost.row_assemblies += 1;
    }

    /// **Where the twenty per-locus buffers' bytes are and how many of them there are** —
    /// the cheap half of the loop's zero-allocation invariant
    /// (`doc/devel/ng/spec/calling_em_loop.md` §13's test 7).
    ///
    /// A `Vec` that **grew** during the loop moves its bytes, so its pointer changes; one
    /// refilled in place does not.
    ///
    /// **⚠ Two of them exchange pointers with every pass, and that is by design**:
    /// [`advance_cohort_expected_copies`](Self::advance_cohort_expected_copies) `mem::swap`s the
    /// cohort's expected copies with the previous pass's rather than copying either. So this
    /// list's *order* depends on the parity of the pass count, and a caller comparing two runs
    /// must compare it as a set — sorted — or compare only runs whose pass counts have the same
    /// parity. Measured: a 2-pass run and a 7-pass run over identical evidence differ in exactly
    /// those two entries, both still `(pointer, length)` pairs the other list holds.
    ///
    /// **What this cannot see is a buffer that was reallocated and came back to the same
    /// address**, which is what a freed block of the same size usually gets — so a `Vec` cloned
    /// once a row leaves this list unchanged after an even number of clones. It cannot see a
    /// temporary allocated and dropped inside a pass either, which leaves no trace anywhere.
    /// **Both halves are counted for real**, by `tests/ng_calling_loop_allocation.rs`, which
    /// installs `dhat`'s counting allocator and reads `total_blocks` across runs at different
    /// pass counts. Measured: reallocating the contaminant frequency table once a row once a
    /// pass leaves every fingerprint here identical and takes that counter from 8 blocks to 24.
    /// It lives in a test binary of its own because a global allocator counts the whole process,
    /// and the lib suite runs its tests in parallel.
    ///
    /// **The row scratches are deliberately not here**, and their absence is the invariant
    /// rather than a gap: each `GenericRowScratch` sizes itself per *sample*, inside the
    /// emission build, so it legitimately grows within a locus when a wider sample arrives. The
    /// emission build happens once, outside the frequency loop, so a pass still allocates
    /// nothing — which is what §13's test 7 claims. Fingerprinting them would make this test
    /// fail on correct code. **The three contaminant buffers above are here**, because they are
    /// refilled inside the loop and sized outside it, which is exactly the property this
    /// measures.
    ///
    /// Test-only, because the pointers are an implementation detail that no run should read.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn buffer_fingerprints(&self) -> Vec<(usize, usize)> {
        fn of<T>(buffer: &[T]) -> (usize, usize) {
            (buffer.as_ptr() as usize, buffer.len())
        }
        vec![
            of(&self.genotype_likelihoods),
            of(&self.prior_row),
            of(&self.posterior_row),
            of(&self.seed_concentration),
            of(&self.sample_concentration),
            of(&self.prior_per_allele_workspace),
            of(&self.cohort_expected_copies),
            of(&self.previous_cohort_expected_copies),
            of(&self.per_sample_expected_copies),
            of(&self.copy_count_log_likelihoods),
            of(&self.allele_count_distribution),
            of(&self.allele_count_distribution_next),
            of(&self.log_allele_count_distribution),
            of(&self.pooled_allele_reads),
            of(&self.wrong_read_charges),
            of(&self.run_sample_of_each_row),
            of(&self.inbreeding_coefficient_by_row),
            of(&self.batch_allele_copies),
            of(&self.contaminant_allele_frequencies),
            of(&self.own_contaminant_row),
            of(&self.expected_copies_by_run_sample),
        ]
    }

    /// **Size this locus's two contaminant tables** — one row per sequencing batch, and the
    /// run's whole sample axis for the copies those rows are summed from.
    ///
    /// **Called only where the locus's own mixture reads a per-batch frequency**, which is an
    /// **ordinary site** in a run that fitted a fraction — a repeat tract's third term is the
    /// fit's frozen allele length variation and reads none. And called after
    /// [`prepare_for_locus`](Self::prepare_for_locus), which un-sizes them. That order is what
    /// makes a locus unable to read another's frequencies: the tables come back empty and every
    /// accessor below refuses them.
    ///
    /// **The run's sample axis and not the locus's rows.** The batching is declared over the
    /// run, so a batch whose samples were all ruled uncallable at this locus still has a row in
    /// the table; summing over rows would leave it unwritten, which
    /// [`fill_batch_allele_copies`] refuses by name. Scattering onto the run's axis gives it a
    /// row of zeros — the same thing the M-step already does with an uncallable sample.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch; on a run with no batches or no samples, which is a run whose
    /// batching went missing rather than a run without one; and on a table longer than a
    /// `usize` can index.
    pub(crate) fn prepare_contaminant_tables(&mut self, batch_count: usize, run_samples: usize) {
        self.assert_prepared();
        assert!(
            batch_count > 0,
            "a contaminated run has at least one sequencing batch, the default being one that \
             holds all of it, so a locus prepared for none is a run whose batching went missing"
        );
        // **The precondition the scatter actually has**, which is not `run_samples >=
        // row_count`: the copies are written at `run_sample × alleles`, so what must hold is
        // that every claimed row names a sample the run has. A one-row locus whose row names
        // sample 5 of a two-sample run passes the count comparison and panics later on a slice
        // range that names neither the row nor the sample.
        self.assert_rows_claimed();
        let highest_named = self
            .run_sample_of_each_row
            .iter()
            .copied()
            .max()
            .expect("a prepared locus has claimed at least one row");
        assert!(
            highest_named < run_samples,
            "this locus's rows name sample {highest_named} and the run has {run_samples}; a row \
             is one of the run's samples, so a row naming one the run does not have is not this \
             run's locus"
        );
        let batch_table_len = batch_count
            .checked_mul(self.allele_count)
            .unwrap_or_else(|| {
                panic!(
                    "a run of {batch_count} sequencing batches over {} alleles needs a contaminant \
                 table longer than a usize can index",
                    self.allele_count
                )
            });
        let run_table_len = run_samples
            .checked_mul(self.allele_count)
            .unwrap_or_else(|| {
                panic!(
                    "a run of {run_samples} samples over {} alleles needs a copies table longer \
                 than a usize can index",
                    self.allele_count
                )
            });
        // **All three take the unwritten sentinel here**, because none of them is read before
        // something writes it: `batch_copy_buffers_mut` zeroes the run-keyed copies in the same
        // call that scatters into them, and the two batch tables are written whole. The zero
        // that matters is that one, and its reason lives beside it.
        resize_and_fill(
            &mut self.expected_copies_by_run_sample,
            run_table_len,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.batch_allele_copies,
            batch_table_len,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.contaminant_allele_frequencies,
            batch_table_len,
            UNWRITTEN_SCRATCH_VALUE,
        );
        resize_and_fill(
            &mut self.own_contaminant_row,
            self.allele_count,
            UNWRITTEN_SCRATCH_VALUE,
        );
        self.batch_count = batch_count;
        self.run_sample_count = run_samples;
    }

    /// **Scatter this pass's per-row expected copies back onto the run's sample axis**, then
    /// hand back that axis and the batch table to sum it into.
    ///
    /// The scatter is here rather than at the call site because it is the join the two axes
    /// meet at: the row map is this scratch's and so is the run-keyed buffer, and a caller
    /// writing the loop itself would be writing this scratch's row arithmetic a second time.
    ///
    /// # Panics
    ///
    /// On a scratch whose contaminant tables were not prepared.
    pub(crate) fn batch_copy_buffers_mut(&mut self) -> BatchCopyBuffers<'_> {
        self.assert_contaminant_tables_prepared();
        self.assert_rows_claimed();
        // **Zeroed here rather than at sizing, and the zero is an answer rather than a slot
        // nobody filled.** Every sample this locus has no row for keeps it, and zero copies is
        // exactly what an uncallable sample contributes to the cohort's expected copies too — so
        // a batch made entirely of samples this locus set aside gets a row of zeros and takes
        // the frequency's no-evidence fallback, rather than a row nobody wrote.
        self.expected_copies_by_run_sample.fill(0.0);
        let allele_count = self.allele_count;
        for (row, &run_sample) in self.run_sample_of_each_row.iter().enumerate() {
            let from = row * allele_count..(row + 1) * allele_count;
            let onto = run_sample * allele_count..(run_sample + 1) * allele_count;
            self.expected_copies_by_run_sample[onto]
                .copy_from_slice(&self.per_sample_expected_copies[from]);
        }
        BatchCopyBuffers {
            expected_copies_by_run_sample: &self.expected_copies_by_run_sample,
            allele_count,
            batch_allele_copies: &mut self.batch_allele_copies,
        }
    }

    /// The three buffers one **sample's** contaminant frequencies are built between.
    ///
    /// # Panics
    ///
    /// On a scratch whose contaminant tables were not prepared, and on a row past what this
    /// locus was prepared for.
    pub(crate) fn contaminant_frequency_buffers_mut(
        &mut self,
        row: usize,
    ) -> ContaminantFrequencyBuffers<'_> {
        self.assert_contaminant_tables_prepared();
        let range = self.allele_row_range(row);
        ContaminantFrequencyBuffers {
            batch_allele_copies: &self.batch_allele_copies,
            own_expected_copies: &self.per_sample_expected_copies[range],
            allele_count: self.allele_count,
            own_contaminant_row: &mut self.own_contaminant_row,
            own_contaminant_batch: &mut self.own_contaminant_batch,
        }
    }

    /// **Fill the frequencies every sample of this pass shares, and lend them out for the pass**
    /// ([`fill_shared_contaminant_allele_frequencies`](crate::calling::likelihood::fill_shared_contaminant_allele_frequencies)).
    /// They leave the scratch so that each row's buffers can be borrowed beside them; give them
    /// back with [`Self::return_shared_contaminant_frequencies`]. Returns them and how many
    /// batches had nothing to read a frequency off.
    ///
    /// # Panics
    ///
    /// On a scratch whose contaminant tables were not prepared.
    pub(crate) fn lend_shared_contaminant_frequencies(&mut self) -> (Vec<f64>, usize) {
        self.assert_contaminant_tables_prepared();
        let mut shared = std::mem::take(&mut self.contaminant_allele_frequencies);
        let no_evidence = crate::calling::likelihood::fill_shared_contaminant_allele_frequencies(
            &self.batch_allele_copies,
            self.allele_count,
            &mut shared,
        );
        (shared, no_evidence)
    }

    /// Take back what [`Self::lend_shared_contaminant_frequencies`] lent.
    pub(crate) fn return_shared_contaminant_frequencies(&mut self, shared: Vec<f64>) {
        self.contaminant_allele_frequencies = shared;
    }

    /// How many sequencing batches this locus's contaminant tables hold — **zero where they
    /// were never prepared**, which is every locus of an uncontaminated run.
    ///
    /// Test-only: a run reads the tables, never their shape.
    #[cfg(test)]
    #[inline]
    #[must_use]
    pub(crate) fn contaminant_batch_count(&self) -> usize {
        self.batch_count
    }

    /// The contaminant allele frequencies **the last row was scored against**, batch-major: the
    /// shared table with that row's own batch row in place of the shared one.
    ///
    /// Test-only: in a run they are read through
    /// [`generic_row_buffers_mut`](Self::generic_row_buffers_mut), one row at a time and while
    /// that row is being assembled, so nothing outside a test ever sees them afterwards.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn contaminant_allele_frequencies(&self) -> Vec<f64> {
        let mut table = self.contaminant_allele_frequencies.clone();
        if !table.is_empty() {
            let from = self.own_contaminant_batch * self.allele_count;
            table[from..from + self.allele_count].copy_from_slice(&self.own_contaminant_row);
        }
        table
    }

    /// Each sequencing batch's expected allele copies at this locus, batch-major, as the last
    /// fill left them. Test-only, for [`contaminant_allele_frequencies`](Self::contaminant_allele_frequencies)'s
    /// reason.
    #[cfg(test)]
    #[must_use]
    pub(crate) fn batch_allele_copies(&self) -> &[f64] {
        &self.batch_allele_copies
    }

    /// Refuse a scratch whose contaminant tables were never sized for this locus.
    ///
    /// **It says which table went missing, which is the whole of what it adds.** The fills do
    /// refuse an empty table on their own — `fill_batch_allele_copies` compares the copies
    /// against `samples × alleles` before anything else — but they refuse it as a shape
    /// mismatch between two slices, at a locus, with nothing in the message about the one call
    /// that was skipped. This names it.
    #[inline]
    fn assert_contaminant_tables_prepared(&self) {
        self.assert_prepared();
        assert!(
            self.batch_count > 0,
            "this locus's contaminant tables were not prepared: `prepare_contaminant_tables` \
             runs after `prepare_for_locus` on a run that fitted a contamination fraction, and \
             nothing else sizes them"
        );
    }

    /// One scratch row's SNP/indel emission cache, which owns its own sizing.
    ///
    /// # Panics
    ///
    /// On a row past what [`prepare_for_locus`](Self::prepare_for_locus) sized this scratch
    /// for.
    #[inline]
    pub fn generic_row_mut(&mut self, row: usize) -> &mut GenericRowScratch {
        assert!(
            row < self.row_count,
            "row {row} is past the {} this scratch was prepared for",
            self.row_count
        );
        &mut self.generic_rows[row]
    }

    #[cfg(test)]
    /// The repeat-tract row's own scratch, which owns its own sizing.
    #[inline]
    pub fn ssr_row_mut(&mut self) -> &mut SsrRowScratch<SsrEmissionScratch> {
        &mut self.ssr_row
    }

    /// **The tract's fitted scoring parameters, to gather into once per repeat tract.**
    ///
    /// Cleared and refilled by
    /// [`TractScoringFits::gather_for_locus`](super::calling::inference::repeat_tract_parameters::TractScoringFits::gather_for_locus),
    /// which is the only thing that writes it.
    #[inline]
    pub(crate) fn tract_fits_mut(
        &mut self,
    ) -> &mut crate::calling::inference::repeat_tract_parameters::TractScoringFits {
        &mut self.tract_fits
    }

    /// The same fits to read — **what the locus's warrant is folded from**, once the rows have
    /// been scored.
    ///
    /// **Separate from [`tract_locus_buffers_mut`](Self::tract_locus_buffers_mut) because the
    /// warrant is read after the scoring rather than during it**, and by then nothing else of
    /// this scratch is borrowed.
    #[inline]
    pub(crate) fn tract_fits(
        &self,
    ) -> &crate::calling::inference::repeat_tract_parameters::TractScoringFits {
        &self.tract_fits
    }

    #[cfg(test)]
    /// How many per-locus slippage re-fit rounds the last locus adopted — the loop's probe
    /// for its own stopping rule, written by the driver after the rounds stop. Zero at the
    /// shipped configuration and at every SNP/indel locus.
    #[inline]
    #[must_use]
    pub fn slippage_refit_rounds(&self) -> u32 {
        self.slippage_refit_rounds
    }

    /// Record how many slippage re-fit rounds this locus adopted — the driver's write behind
    /// [`Self::slippage_refit_rounds`], once per locus after its rounds stop.
    #[inline]
    pub(crate) fn record_slippage_refit_rounds(&mut self, rounds: u32) {
        self.slippage_refit_rounds = rounds;
    }

    /// **Everything a repeat tract's rows are scored from, borrowed from one scratch at once.**
    ///
    /// **It exists because the contexts the row takes borrow the fits, and the row writes into
    /// this same scratch.** The contexts are built once per tract, from `fits`, and every row
    /// then needs them alive while `ssr_row` and the likelihood table are written — three
    /// borrows of three different fields, which no sequence of one-field accessors can hold at
    /// once. **The borrows are disjoint by construction**: every field below is a different
    /// field of this scratch.
    ///
    /// **The whole table rather than one row**, unlike the SNP/indel path's, because the tract's
    /// row scratch is one buffer reused across rows rather than one per row — so the loop over
    /// rows lives inside this borrow instead of outside it.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch, and on a row map that is not one entry per prepared row —
    /// both of which say the locus's shape was never fixed. **The row map is checked here
    /// rather than where it is read**, because the walk over rows reads it inside this borrow,
    /// where no accessor of this scratch can be called.
    #[inline]
    pub(crate) fn tract_locus_buffers_mut(&mut self) -> TractLocusBuffers<'_, SsrEmissionScratch> {
        self.assert_rows_claimed();
        TractLocusBuffers {
            fits: &self.tract_fits,
            run_sample_of_each_row: &self.run_sample_of_each_row,
            row_scratch: &mut self.ssr_row,
            genotype_likelihoods: &mut self.genotype_likelihoods,
            genotype_count: self.genotype_count,
        }
    }

    /// Every buffer scoring **one sample** at this locus reads or writes, handed out
    /// together.
    ///
    /// **It exists because four of them have to be live at once and they are fields of one
    /// scratch.** Scoring a sample builds its concentration from three buffers into a
    /// fourth, its log-prior from that concentration and a working buffer into a fifth, its
    /// posterior from that prior and its likelihood row into a sixth, and finally its own
    /// expected copies from the posterior. Reaching for those one accessor at a time does
    /// not compile — each borrows the whole scratch — and the way round it that needs no
    /// new type is to copy a buffer, which would put an allocation back into a pass whose
    /// whole shape exists to have none.
    ///
    /// **The borrows are disjoint by construction**: every borrowed field below is a
    /// different field of this scratch, so nothing here can alias.
    ///
    /// **`_mut` because it is the mutable half of five accessor pairs at once**, not the
    /// shared half of a sixth — five of the eight borrows it hands out are `&mut`, and the
    /// suffix is the only thing on the name that says a caller needs a mutable scratch.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch, or a sample past the count
    /// [`prepare_for_locus`](Self::prepare_for_locus) was given.
    #[inline]
    #[must_use]
    pub(crate) fn sample_scoring_buffers_mut(&mut self, sample: usize) -> SampleScoringBuffers<'_> {
        let genotype_row = self.genotype_row_range(sample);
        let allele_row = self.allele_row_range(sample);
        SampleScoringBuffers {
            sample,
            seed_concentration: &self.seed_concentration,
            cohort_expected_copies: &self.cohort_expected_copies,
            genotype_likelihoods: &self.genotype_likelihoods[genotype_row],
            sample_concentration: &mut self.sample_concentration,
            prior_per_allele_workspace: &mut self.prior_per_allele_workspace,
            prior_row: &mut self.prior_row,
            posterior_row: &mut self.posterior_row,
            sample_expected_copies: &mut self.per_sample_expected_copies[allele_row],
        }
    }

    /// The two buffers the **M-step** reads and writes: every sample's own expected allele
    /// copies, and the cohort row their sum goes into.
    ///
    /// **A second bundle for the same reason as the first** — the per-sample copies and the
    /// cohort's are two fields of one scratch, and reaching for them one accessor at a time
    /// does not compile. Measured while B1 was in review: an M-step written against the
    /// per-buffer accessors gives two `error[E0502]`.
    ///
    /// **The whole per-sample table, not one row.** The M-step's defining property is that it
    /// adds the samples up *in the run's fixed order*, so it takes the table and walks it
    /// itself rather than being handed a row at a time by a caller whose loop could reorder.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub(crate) fn cohort_summing_buffers_mut(&mut self) -> CohortSummingBuffers<'_> {
        self.assert_prepared();
        CohortSummingBuffers {
            sample_count: self.row_count,
            per_sample_expected_copies: &self.per_sample_expected_copies,
            cohort_expected_copies: &mut self.cohort_expected_copies,
        }
    }

    /// The artifact summary's per-allele read totals, to zero and refill.
    ///
    /// **Handed out mutably and never read back through a shared accessor**, because the
    /// one caller that fills it is the one caller that reads it, in the same call: the
    /// final pass pools every callable sample's reads onto the alleles, calls the samples,
    /// picks the primary alternative from the totals and the calls, and is done with them
    /// (`doc/devel/ng/spec/calling_quality.md` §6.3). There is nothing for a later stage to
    /// come back for.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub(crate) fn pooled_allele_reads_mut(&mut self) -> &mut [u64] {
        self.assert_prepared();
        &mut self.pooled_allele_reads
    }

    /// The five buffers the **site quality** reads and writes: the likelihood table the loop
    /// built, and the four the fold needs.
    ///
    /// **A third bundle, for the third reason of the same kind** — five fields of one
    /// scratch, and the fold swaps two of them, which needs both to be borrowed at once from
    /// the same call.
    ///
    /// **Called once per locus, after the loop has stopped.** The likelihood table it hands
    /// over is the one the loop leaves behind, which is what makes the site quality cost no
    /// emission evaluations of its own (`doc/devel/ng/spec/calling_quality.md` §3.2) — and
    /// where a fraction was fitted, the one assembled against the frequencies the loop settled
    /// on, so the quality and the genotypes beside it come from the same table.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch.
    #[inline]
    #[must_use]
    pub fn site_quality_buffers_mut(&mut self) -> SiteQualityBuffers<'_> {
        self.assert_prepared();
        SiteQualityBuffers {
            sample_count: self.row_count,
            genotype_likelihoods: &self.genotype_likelihoods,
            copy_count_log_likelihoods: &mut self.copy_count_log_likelihoods,
            allele_count_distribution: &mut self.allele_count_distribution,
            allele_count_distribution_next: &mut self.allele_count_distribution_next,
            log_allele_count_distribution: &mut self.log_allele_count_distribution,
            count_prior: &mut self.count_prior,
        }
    }

    /// Refuse a scratch that was never sized for a locus.
    ///
    /// # Panics
    ///
    /// If [`prepare_for_locus`](Self::prepare_for_locus) has not been called. Every buffer
    /// of a freshly-made scratch is empty, so without this a pass folding into
    /// [`cohort_expected_copies_mut`](Self::cohort_expected_copies_mut) would run zero
    /// iterations, write nothing, and leave the cohort's copies summing to a plausible
    /// `0.0` — which the fill value cannot catch, because an unprepared buffer has no slots
    /// to fill.
    #[inline]
    fn assert_prepared(&self) {
        assert!(
            self.row_count > 0,
            "this scratch has not been prepared for a locus: call prepare_for_locus first"
        );
    }

    /// Where one sample's row of the genotype-likelihood table starts and ends.
    #[inline]
    fn genotype_row_range(&self, sample: usize) -> std::ops::Range<usize> {
        self.row_range(sample, self.genotype_count, "genotype-likelihood")
    }

    /// Where one sample's row of the per-sample expected copies starts and ends.
    #[inline]
    fn allele_row_range(&self, sample: usize) -> std::ops::Range<usize> {
        self.row_range(sample, self.allele_count, "per-sample expected copies")
    }

    /// Where one sample's row of a flat sample-major table starts and ends.
    ///
    /// **Private, and reached only through the two typed wrappers**, so no caller picks the
    /// row width by hand — the two widths are the genotype count and the allele count, and
    /// at a diploid biallelic locus they are 3 and 2, so a wrong one is a legal range over
    /// the wrong table.
    ///
    /// # Panics
    ///
    /// On an unprepared scratch, or a sample past the count
    /// [`prepare_for_locus`](Self::prepare_for_locus) was given — a caller that walked one
    /// sample too far would otherwise read the next sample's row, or, at the last sample,
    /// panic somewhere that says nothing about which sample it was.
    #[inline]
    fn row_range(&self, sample: usize, width: usize, table: &str) -> std::ops::Range<usize> {
        self.assert_prepared();
        assert!(
            sample < self.row_count,
            "sample {sample} is past the {} this scratch was prepared for, indexing the \
             {table} table",
            self.row_count
        );
        let start = sample * width;
        start..start + width
    }
}

/// Every buffer one sample's turn through a pass touches, borrowed from one
/// [`CallingScratch`] at once.
///
/// **Every run-time construction goes through
/// [`CallingScratch::sample_scoring_buffers_mut`]**, which is where the flat sample-major
/// tables are sliced, so the row arithmetic has one spelling. See that method for why the
/// buffers travel together rather than one at a time.
///
/// **The fields are public because the tests build one by hand**, deliberately: the only way
/// to hand the scorer a row of the wrong width — the mis-shape its release-held checks exist
/// for — is to reach past the accessor, which cannot produce one. Nothing in a run should.
/// The shapes are not this type's to enforce: they are checked where they are used, in
/// `score_one_sample`, `fill_sample_concentration` and `PriorRow::new`.
///
/// Of the eight buffers, three are read-only — what the locus's setup and the *previous* pass
/// left — and four are this sample's share of what this pass produces. The eighth is both:
/// this sample's expected copies are read as the leave-one-out term the prior needs, and then
/// overwritten with what this pass's posterior implies. The `sample` index beside them is not
/// a buffer; it is there so a panic can say which sample was being scored.
#[derive(Debug)]
pub(crate) struct SampleScoringBuffers<'a> {
    /// Which sample of the run these buffers belong to — carried so a panic raised while
    /// scoring can name it. At a thousand samples that is the first thing the reader of a
    /// message wants, and no other field records it.
    pub sample: usize,
    /// The locus's seed concentration — one entry per allele, the same for every sample,
    /// and what the prior falls back to where the cohort says nothing.
    pub seed_concentration: &'a [f64],
    /// The cohort's expected allele copies as the previous pass left them, this sample's
    /// own contribution included.
    pub cohort_expected_copies: &'a [f64],
    /// This sample's row of the genotype-likelihood table — one entry per candidate
    /// genotype, assembled from emissions the locus computed once.
    pub genotype_likelihoods: &'a [LogProb],
    /// One entry per allele, to fill: this sample's own concentration — the seed plus what
    /// the **other** samples showed here.
    pub sample_concentration: &'a mut [f64],
    /// One entry per allele of working space, to hand the genotype prior. Its contents on
    /// entry are ignored and on exit unspecified.
    pub prior_per_allele_workspace: &'a mut [f64],
    /// One entry per candidate genotype, to fill: this sample's log-prior.
    pub prior_row: &'a mut [LogProb],
    /// One entry per candidate genotype, to fill: this sample's posterior over genotypes.
    pub posterior_row: &'a mut [f64],
    /// One entry per allele: this sample's own expected allele copies — **read first** as
    /// the leave-one-out term, then **overwritten** with this pass's answer.
    pub sample_expected_copies: &'a mut [f64],
}

/// The three buffers the **SNP/indel row builder** reads and writes, borrowed from one
/// [`CallingScratch`] in one call.
///
/// **A bundle for the reason the other three are**: they are three fields of one scratch, and
/// reaching for them an accessor at a time does not compile. The error spreads are shared by
/// every sample and filled once per locus; the other two are this row's.
#[derive(Debug)]
pub(crate) struct GenericRowBuffers<'a> {
    /// The locus's per-allele wrong-read charge table, filled once before the first row.
    pub(crate) wrong_read_charges: &'a [WrongReadCharge],
    /// This row's own batch's contaminant frequencies, with the sample left out — the one row
    /// of its table that is its own; the rest are the locus's shared table. Empty on a run the
    /// fit found no contamination in, where the mixture is absent.
    pub(crate) own_contaminant_row: &'a [f64],
    /// This row's own emission cache, which it sizes itself.
    pub(crate) row_scratch: &'a mut GenericRowScratch,
    /// This sample's row of the genotype-likelihood table, to fill.
    pub(crate) genotype_likelihoods: &'a mut [LogProb],
}

/// **Everything one repeat tract's rows are scored from**, borrowed from one
/// [`CallingScratch`] at once — [`CallingScratch::tract_locus_buffers_mut`] is the only way to
/// make one.
///
/// **It hands out the whole likelihood table rather than one row**, which is the difference from
/// [`GenericRowBuffers`]: a tract's row scratch is one reused buffer rather than one per row, so
/// the walk over rows happens inside this borrow.
#[derive(Debug)]
pub(crate) struct TractLocusBuffers<'a, SsrEmissionScratch> {
    /// This tract's fitted scoring parameters — what the contexts are built from, and what the
    /// locus's warrant is folded from.
    pub(crate) fits: &'a crate::calling::inference::repeat_tract_parameters::TractScoringFits,
    /// Which sample of the run each scratch row holds — **the map the walk over rows needs and
    /// cannot ask for separately**, since asking would be a second borrow of this scratch while
    /// the three below are held.
    ///
    /// **At a repeat tract this is provably the run's sample order unchanged**, because a tract
    /// rules no sample out — so no test of the tract path exercises the difference between a
    /// row and the sample it holds. It is carried rather than assumed so that the day a tract
    /// can set a sample aside, the walk that reads it needs no change.
    pub(crate) run_sample_of_each_row: &'a [usize],
    /// The repeat-tract row's working memory, reused across rows.
    pub(crate) row_scratch: &'a mut SsrRowScratch<SsrEmissionScratch>,
    /// `rows × genotypes`, row-major: the whole genotype-likelihood table, to fill.
    pub(crate) genotype_likelihoods: &'a mut [LogProb],
    /// How many candidate genotypes this locus has — the table's stride.
    pub(crate) genotype_count: usize,
}

/// The buffers each **sequencing batch's** expected allele copies are summed between, borrowed
/// from one [`CallingScratch`] at once — [`CallingScratch::batch_copy_buffers_mut`] is the only
/// way to make one, because it is what performs the scatter the first field's name implies.
#[derive(Debug)]
pub(crate) struct BatchCopyBuffers<'a> {
    /// This pass's expected allele copies on the **run's** sample axis, sample-major, zero at
    /// every sample this locus has no row for.
    pub(crate) expected_copies_by_run_sample: &'a [f64],
    /// How many alleles this locus is called over — the stride of both tables.
    pub(crate) allele_count: usize,
    /// One row per sequencing batch, to fill.
    pub(crate) batch_allele_copies: &'a mut [f64],
}

/// The buffers **one sample's** contaminant allele frequencies are built between, borrowed from
/// one [`CallingScratch`] at once.
#[derive(Debug)]
pub(crate) struct ContaminantFrequencyBuffers<'a> {
    /// Every batch's expected allele copies at this locus, batch-major — this sample's own
    /// among them, which is what the fill takes back out.
    pub(crate) batch_allele_copies: &'a [f64],
    /// This sample's own expected allele copies.
    pub(crate) own_expected_copies: &'a [f64],
    /// How many alleles this locus is called over — the stride of both tables.
    pub(crate) allele_count: usize,
    /// This sample's own batch row, to fill: that batch's frequencies with the sample left out.
    pub(crate) own_contaminant_row: &'a mut [f64],
    /// Which batch the row is, recorded for [`CallingScratch::contaminant_allele_frequencies`].
    pub(crate) own_contaminant_batch: &'a mut usize,
}

/// The two buffers the **M-step** works between, borrowed from one [`CallingScratch`] at once.
///
/// **The M-step** is the half of one pass that turns every sample's genotype probabilities back
/// into the cohort-wide summary the next pass conditions on: it adds up how many copies of each
/// allele the samples are expected to carry. Its counterpart is the E-step, which goes the other
/// way — [`SampleScoringBuffers`] is that half's bundle
/// (`doc/devel/ng/spec/calling_em_loop.md` §2).
///
/// **Made only by [`CallingScratch::cohort_summing_buffers_mut`]** in a run; the tests build one
/// by hand to hand the sum a mis-shaped table, which the accessor cannot produce.
#[derive(Debug)]
pub(crate) struct CohortSummingBuffers<'a> {
    /// How many samples the table below holds — the run's whole sample count, which is what
    /// the scratch was prepared for.
    ///
    /// **Deliberately redundant against the table's own length, because that redundancy is
    /// what makes the shape check a check.** Derive the count from the length instead and any
    /// table divides evenly by the allele count, so a two-row table presented as three samples
    /// is accepted and the cohort comes back short by a sample with nothing raised. Measured
    /// during B2's review: with the count derived, the suite stays green and the
    /// wrong-size test starts accepting that table and returning `[2.0, 2.0]`.
    pub sample_count: usize,
    /// `samples × alleles`, sample-major: each sample's own expected allele copies, as the
    /// E-step of this pass left them. Read whole, so that the order of the sum is this
    /// function's and not its caller's.
    pub per_sample_expected_copies: &'a [f64],
    /// One entry per allele, to fill: the cohort's expected allele copies, the sum of the
    /// rows above.
    pub cohort_expected_copies: &'a mut [f64],
}

/// Resize a scratch buffer and overwrite every entry, including the ones that were already
/// there.
///
/// **`Vec::resize` alone is the bug this exists to avoid**: it leaves the leading entries
/// untouched, so a locus with the same shape as the last one would silently reuse the last
/// one's values — a wrong prior or a wrong likelihood at every sample, with nothing failing.
fn resize_and_fill<T: Copy>(buffer: &mut Vec<T>, len: usize, fill: T) {
    buffer.clear();
    buffer.resize(len, fill);
}

/// One sample's call at one locus.
///
/// **Two outcomes, and the second is not a low-confidence version of the first.** A sample
/// the candidate step declared uncallable has no genotype at all: the locus is called over a
/// set of alleles that cannot represent what the sample carries, so it is set aside before
/// the first pass, scored against nothing, and emitted as missing
/// (`doc/devel/ng/spec/calling_em_loop.md` §5.0, §9). Emission must not write that as a
/// genotype with a poor quality beside it, and an enum is what stops it: there is no
/// quality to read.
///
/// **A sample with no reads at the locus is [`Self::Called`], not [`Self::Missing`]**, and
/// the distinction is the one a reader is most likely to get backwards. Such a sample scores
/// every genotype alike, so the prior decides it alone and a genotype comes out — which is
/// the right answer rather than a special case
/// (`doc/devel/ng/spec/calling_em_loop.md` §7). [`Self::Missing`] is not *missing data*; it
/// is *the caller declined to invent a genotype over a set that cannot hold this sample's
/// allele*.
///
/// **Why the type changed shape here.** `doc/devel/ng/arch/calling_em_loop.md` §2 sketches
/// this as a struct of a genotype and a quality, which cannot express the missing case, and
/// [`Genotype::new`] refuses an empty multiset by design — *an empty multiset is not a
/// haploid call, it is a sample with no genome*. The spec's §5.0 records that the ruling had
/// a producer and no carrier and that the loop's own plan adds one; this is it.
#[derive(Clone, PartialEq, Debug)]
pub enum SampleGenotypeCall {
    /// The caller reached a genotype for this sample.
    Called {
        /// Which alleles this sample carries, one per copy of its genome.
        genotype: Genotype,
        /// How sure the caller is of that genotype — the posterior-derived `GQ` of the
        /// loop's last pass, how much of the sample's genotype probability the winning
        /// genotype took. Step 13's quality model **refines** this number; it does not
        /// replace it (`doc/devel/ng/arch/calling_em_loop.md` §2).
        genotype_quality: Phred,
        /// **Whether this sample's own reads said nothing about which genotype it has** —
        /// its genotype likelihoods were flat, so every genotype was equally probable
        /// under them and the prior decided alone.
        ///
        /// **It is what turns a called sample into a `./.`** at emission
        /// (`doc/devel/ng/spec/vcf_output.md` §7.1), and it is recorded **here**, by
        /// whoever scored the sample, because nothing downstream can recover it: the
        /// likelihoods live in per-sample scratch the loop overwrites at the next locus.
        /// `vcf::assemble`'s own module note named this as the one input it could not be
        /// given; this is the answer to it.
        ///
        /// **It must be the likelihood and not the posterior.** A sample with no reads is
        /// scored by the loop and comes back with a genotype, because the prior decides it
        /// alone — and where the fitted frequency is low that posterior is sharply peaked,
        /// so no threshold on the confidence would catch this sample. The likelihood is
        /// what the reads said; the posterior is what the reads said plus what the cohort
        /// assumed.
        ///
        /// **True is common and is not a defect.** On six tomato accessions over 400 kb,
        /// every locus is called on every sample of the run and most samples cover most
        /// loci — but a sample that covered nothing at a locus is scored there all the
        /// same, and this is what says so.
        reads_were_uninformative: bool,
        /// The sample's likelihood and posterior rows over every genotype of the table the locus
        /// was called with — the `PL` and `GP` the VCF writes (`vcf_output.md` §7.2). Taken in
        /// the final pass, where both exist for one sample at a time.
        genotype_rows: GenotypeRows,
    },
    /// The candidate step ruled this sample uncallable at this locus, so it took no part in
    /// the loop and emission writes its `GT` as missing.
    ///
    /// **This happens on the SNP/indel path only.** A repeat tract sets no sample aside
    /// (`doc/devel/ng/spec/calling_em_loop.md` §5.0.1), and [`LocusInference::new`] refuses
    /// a locus that carries one there.
    Missing,
}

#[cfg(test)]
impl SampleGenotypeCall {
    /// **Test fixtures only**: a called sample's rows replaced by flat ones — every genotype equally
    /// likely and equally probable — over `allele_count` alleles at the genotype's ploidy, so a
    /// fixture built without a calling pass still carries rows the VCF assembly can index.
    #[must_use]
    pub fn with_flat_rows(self, allele_count: usize) -> Self {
        match self {
            Self::Called {
                genotype,
                genotype_quality,
                reads_were_uninformative,
                genotype_rows: _,
            } => {
                let ploidy = Ploidy::try_new(
                    u8::try_from(genotype.alleles().len()).expect("a fixture's ploidy fits a u8"),
                )
                .expect("a fixture's genotype names at least one copy");
                let genotypes = GenotypeTable::build(ploidy, allele_count).genotype_count();
                Self::Called {
                    genotype,
                    genotype_quality,
                    reads_were_uninformative,
                    genotype_rows: GenotypeRows {
                        ploidy,
                        log_likelihoods: vec![0.0; genotypes].into(),
                        posteriors: vec![1.0 / genotypes as f64; genotypes].into(),
                    },
                }
            }
            Self::Missing => Self::Missing,
        }
    }
}

/// **Test fixtures only**: rows with no entries, for a call that is never assembled into a record.
#[cfg(test)]
#[must_use]
pub fn placeholder_genotype_rows() -> GenotypeRows {
    GenotypeRows {
        ploidy: Ploidy::try_new(2).expect("two copies is a ploidy"),
        log_likelihoods: Box::new([]),
        posteriors: Box::new([]),
    }
}

impl SampleGenotypeCall {
    /// The alleles called, or `None` where the sample was set aside.
    #[inline]
    #[must_use]
    pub fn genotype(&self) -> Option<&Genotype> {
        match self {
            Self::Called { genotype, .. } => Some(genotype),
            Self::Missing => None,
        }
    }

    #[cfg(test)]
    /// How sure the caller is, or `None` where the sample was set aside — which carries no
    /// quality, because no genotype was scored to have one.
    #[inline]
    #[must_use]
    pub fn score_best_genotype(&self) -> Option<Phred> {
        match self {
            Self::Called {
                genotype_quality, ..
            } => Some(*genotype_quality),
            Self::Missing => None,
        }
    }

    /// Whether this sample's `GT` is written as missing.
    #[inline]
    #[must_use]
    pub fn is_missing(&self) -> bool {
        matches!(self, Self::Missing)
    }
}

/// **What a repeat tract's call rested on** — the four things about a tract's parameters that
/// the genotypes cannot be asked.
///
/// A tract is scored under two fitted numbers per `(read group, candidate)` pair — how often
/// the copying steps before sequencing add or drop a whole repeat, and how often a base is
/// misread — and under a prior shape saying which lengths are plausible. Each of the three can
/// fall back to a stated constant, and **a call resting on measurements and one resting on
/// constants are different claims about the same reads**
/// (`doc/devel/ng/spec/population_diversity.md` §1). Nothing in the called genotype says which
/// it was, so it travels here.
///
/// # The counts are per `(read group, candidate)` cell, and the table covers the whole run
///
/// A tract's parameter table is `read groups × candidates`, over **every** read group of the
/// run rather than the ones whose reads reached this tract — the read likelihood's own table is
/// indexed by read-group identifier directly. So on a run of many libraries a tract can report
/// defaulted cells on account of a library that sent it nothing. That is the conservative
/// direction and it is what [`Self::scoring_cells`] counts against.
///
/// # Why the fields are read through accessors
///
/// Three of the six are counts of the same type, and two of them are *shares of* another: a
/// record claiming seven defaulted cells out of six, or more cells defaulted by an unknown
/// library than defaulted at all, is arithmetic nonsense that would print as an answer.
/// [`Self::new`] refuses those, and public fields would let a caller build one around it.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct RepeatTractProvenance {
    /// **Which rung of the tract ladder this tract's prior shape came from** — the stratum's own
    /// fitted allele length variation, its own shares at the stated concentration, its motif
    /// period's pooled tracts, or a flat shape at the stated concentration
    /// (`doc/devel/ng/spec/population_diversity.md` §4.4).
    ///
    /// Every rung answers, so this is never *no shape*; what it says is how well founded the
    /// shape was. The SNP/indel path has a ladder of its own with different rungs, resolved once
    /// per run rather than per locus, and none of its rungs can appear here — which is why this
    /// whole value is absent there rather than this field being optional.
    ///
    /// *(It replaced a `seed_diversity_unreachable` flag, whose whole subject was a failure that
    /// no longer exists: the prior used to scale a constructed geometric shape to reproduce a
    /// measured gene diversity, which was impossible above a ceiling the shape itself set.
    /// Seeding from the fit asserts no such scaling — `population_diversity.md` §4.2.)*
    allele_length_variation_rung: AlleleLengthVariationRung,
    /// How many `(read group, candidate)` cells this tract was scored over — the denominator the
    /// three counts below are read against.
    scoring_cells: usize,
    /// How many of them took the shipped stutter row because the slippage fit had no numbers for
    /// them.
    ///
    /// **Ordinary rather than an error**: a candidate several repeats from its tract's reference
    /// length lands in no fitted stratum on perfectly good data.
    cells_with_no_fitted_slippage: usize,
    /// Of those, how many were defaulted because **the slippage fit does not describe this run's
    /// read groups** — a library the pre-pass never saw, or a slippage map naming more libraries
    /// than the fit was run over.
    ///
    /// **This absence means the parameters and the reads came from different runs**, where the
    /// other one is a candidate sitting off the fitted range on perfectly good data. That is why
    /// it is counted apart rather than folded in: anything above zero is a run to look at, not a
    /// tract.
    ///
    /// **It sees the slippage fit only, and the substitution rates can fail the same way
    /// unseen.** Those are looked up in a plain map keyed by `(read group, stratum, ploidy)`
    /// whose absence carries no reason
    /// ([`FrozenParameters::ssr_substitution_rate_at`]), so a rate map fitted over a different
    /// set of libraries lands in the count below, indistinguishable from a stratum that was
    /// simply never fitted. Splitting that count the same way means giving that lookup a typed
    /// absence, which is a change to the parameter side and is recorded as open rather than made
    /// here.
    cells_whose_read_group_the_fit_does_not_describe: usize,
    /// How many cells found no fitted substitution rate and took the stated constant, which is
    /// defined as the SNP/indel path's own default so that a run cannot default its two error
    /// parameters to two different guesses.
    ///
    /// **One count, two causes, and nothing here separates them** — see the field above.
    cells_with_no_fitted_substitution_rate: usize,
    /// **Whether the mixture's third term was built** — how common each reachable tract length
    /// is in the contaminating population.
    ///
    /// True exactly where the run's parameter fit found a contamination fraction, since that is
    /// what decides whether a tract gets the three-term form or the two
    /// (`doc/devel/ng/spec/read_likelihoods.md` §4.5.1). It is a run-wide condition read at the
    /// tract, and it is here rather than only in the run's report because a reader holding one
    /// locus's record should not have to fetch the run's to know whether its reads were shared
    /// out with a contaminant.
    contaminant_term_was_built: bool,
}

impl RepeatTractProvenance {
    /// One repeat tract's record, with its three counts checked against each other.
    ///
    /// # Panics
    ///
    /// Held in release, because each is a record that would print as an answer
    /// (`doc/devel/ng/spec/calling_em_loop.md` §8):
    ///
    /// - **more cells with no fitted slippage than cells.** The count is over the same table its
    ///   denominator measures, so a larger one means the two came from different tracts.
    /// - **more cells defaulted by a library the fit does not describe than defaulted at all.**
    ///   The first is documented as a share of the second and is counted inside its own arm; a
    ///   larger one means the two counters were swapped, which is the one mistake a reader of
    ///   these numbers could not detect.
    /// - **more cells with no fitted substitution rate than cells**, for the first reason again
    ///   on the other parameter.
    #[must_use]
    pub fn new(
        allele_length_variation_rung: AlleleLengthVariationRung,
        scoring_cells: usize,
        cells_with_no_fitted_slippage: usize,
        cells_whose_read_group_the_fit_does_not_describe: usize,
        cells_with_no_fitted_substitution_rate: usize,
        contaminant_term_was_built: bool,
    ) -> Self {
        assert!(
            cells_with_no_fitted_slippage <= scoring_cells,
            "a tract cannot have {cells_with_no_fitted_slippage} of its {scoring_cells} cells \
             fall back on the slippage fit: the count and its denominator are over one table, \
             so a larger count came from a different tract"
        );
        assert!(
            cells_whose_read_group_the_fit_does_not_describe <= cells_with_no_fitted_slippage,
            "{cells_whose_read_group_the_fit_does_not_describe} cells were defaulted by a \
             library the fit does not describe and only {cells_with_no_fitted_slippage} were \
             defaulted at all: the first is a share of the second, so this is the two counters \
             swapped"
        );
        assert!(
            cells_with_no_fitted_substitution_rate <= scoring_cells,
            "a tract cannot have {cells_with_no_fitted_substitution_rate} of its \
             {scoring_cells} cells fall back on the substitution rate: the count and its \
             denominator are over one table"
        );
        Self {
            allele_length_variation_rung,
            scoring_cells,
            cells_with_no_fitted_slippage,
            cells_whose_read_group_the_fit_does_not_describe,
            cells_with_no_fitted_substitution_rate,
            contaminant_term_was_built,
        }
    }

    /// Which rung of the tract ladder this tract's prior shape came from — see
    /// [`Self::allele_length_variation_rung`].
    #[inline]
    #[must_use]
    pub fn allele_length_variation_rung(&self) -> AlleleLengthVariationRung {
        self.allele_length_variation_rung
    }

    /// How many `(read group, candidate)` cells this tract was scored over.
    #[inline]
    #[must_use]
    pub fn scoring_cells(&self) -> usize {
        self.scoring_cells
    }

    /// How many of them took the shipped stutter row — see
    /// [`Self::cells_with_no_fitted_slippage`].
    #[inline]
    #[must_use]
    pub fn cells_with_no_fitted_slippage(&self) -> usize {
        self.cells_with_no_fitted_slippage
    }

    /// Of those, how many because the slippage fit does not describe this run's read groups —
    /// see [`Self::cells_whose_read_group_the_fit_does_not_describe`].
    #[inline]
    #[must_use]
    pub fn cells_whose_read_group_the_fit_does_not_describe(&self) -> usize {
        self.cells_whose_read_group_the_fit_does_not_describe
    }

    /// How many took the stated substitution rate — see
    /// [`Self::cells_with_no_fitted_substitution_rate`].
    #[inline]
    #[must_use]
    pub fn cells_with_no_fitted_substitution_rate(&self) -> usize {
        self.cells_with_no_fitted_substitution_rate
    }

    /// Whether the mixture's third term was built — see [`Self::contaminant_term_was_built`].
    #[inline]
    #[must_use]
    pub fn contaminant_term_was_built(&self) -> bool {
        self.contaminant_term_was_built
    }
}

/// What calling produces at one locus: the alleles it settled on, every sample's call,
/// and the evidence for how it got there.
///
/// **How the answer was reached travels with the answer, because nothing downstream can
/// reconstruct it.** A genotype from a loop that ran out of passes is a different claim
/// from one that settled; a parameter that was defaulted is a guess where a fitted one
/// is a measurement; a prior seed that could not reach its target is a prior that was
/// bent to fit. Each of those is a field here, and each would otherwise be invisible by
/// the time the record is written (`doc/devel/ng/spec/calling_em_loop.md` §6,
/// `doc/devel/ng/arch/calling_priors.md` §5).
///
/// Plain data, and six of its ten fields are public because every consumer reads
/// them. **The four that are not are private for two different reasons.**
/// [`Self::alleles`], [`Self::cohort_expected_copies`] and
/// [`Self::artifact_test_counts`] are all indexed by, or parallel to, the allele table: a
/// public allele table would let a consumer widen it with [`CandidateAlleles::admit`]
/// against unchanged copies, or leave an artifact summary naming an allele the locus no
/// longer holds — breaking, after construction, the pairings [`Self::new`] checks at it.
/// [`Self::site_quality`] is private for a different rule, which its own comment gives:
/// **there is one quality field, and nothing may read it between the worker that writes
/// the baseline and the stage that overwrites it with the corrected value**
/// (`doc/devel/ng/spec/calling_quality.md` §3.5).
#[derive(Clone, PartialEq, Debug)]
pub struct LocusInference {
    /// Where on the reference this locus is.
    pub region: GenomeRegion,
    /// The alleles it was finally called over — after any discovery round admitted
    /// alleles the first pass missed, and after the prune dropped the ones no sample
    /// carried. Read through [`Self::alleles`].
    alleles: CandidateAlleles,
    /// One call per sample, **in the run's sample order** — the same order every
    /// per-sample slice the loop was given is indexed by. That one order for the whole
    /// run is what makes the loop's fixed-order sum over samples reproducible
    /// (`doc/devel/ng/spec/calling_em_loop.md` §8).
    pub per_sample: Vec<SampleGenotypeCall>,
    /// The cohort's expected allele copies, one entry per allele of [`Self::alleles`].
    ///
    /// **Handed on rather than left to be recomputed**: deriving it downstream from the
    /// called genotypes gives a different number, because a call has already thrown away
    /// the uncertainty these counts carry, and site filtering and emission both read it
    /// (`doc/devel/ng/spec/calling_em_loop.md` §9). Read through
    /// [`Self::cohort_expected_copies`].
    cohort_expected_copies: ExpectedAlleleCopies,
    /// Whether the allele frequencies stopped moving, or the loop simply ran out of
    /// passes.
    ///
    /// **`false` is emitted, never dropped and never fatal.** One hard locus must not
    /// kill a cohort run — production retired its non-convergence error for exactly this
    /// reason — but the flag has to reach the output, because a genotype from a loop that
    /// did not settle is a weaker claim than one from a loop that did, and nothing
    /// downstream can tell them apart otherwise
    /// (`doc/devel/ng/spec/calling_em_loop.md` §6).
    pub converged: bool,
    /// How many passes the frequency loop took. At least one.
    ///
    /// **One pass is a legitimate answer, not a defect.** The loop starts with an E-step
    /// that uses no prior at all — the reads alone — and only then begins iterating
    /// (`doc/devel/ng/spec/calling_em_loop.md` §3), so the first pass already has a
    /// previous estimate to be compared against and can settle against it.
    ///
    /// Carried as an instrument, not as a diagnostic: the pass cap and the convergence
    /// threshold are both inherited from production and neither has been measured on
    /// this caller's range — a thousand samples, or three reads a position — so the
    /// distribution of this number across a real panel is what will set them
    /// (`doc/devel/ng/spec/calling_em_loop.md` §12, question 4).
    pub passes: u32,
    /// The weakest warrant among every parameter that entered this locus's arithmetic.
    ///
    /// The scoring contexts propagate it rather than branching on it: each carries the
    /// weakest [`Provenance`] of the parameters that reached it, and the loop copies the
    /// weakest of those onto the locus (`doc/devel/ng/arch/read_likelihoods.md` §1.4). A
    /// consumer that treats a defaulted error rate like a fitted one is the failure this
    /// exists to prevent.
    ///
    /// **Which warrant is weakest is settled**, and this field is the first consumer of the
    /// answer: [`Provenance::weaker_of`] combines two along the ladder
    /// `parameter_estimation`'s own documentation states — fitted here, borrowed from the
    /// sample's other read groups, supplied, defaulted — so a value the run *supplied* ranks
    /// below a borrowed fit, on the ground that a handed-in number says nothing about this
    /// data where a borrowed one is at least a measurement of a neighbouring grain.
    ///
    /// **The two paths fold two different lists, because they read two different
    /// parameters.** A SNP/indel locus folds the calibrations of the read groups whose reads
    /// reached it. A repeat tract folds the stutter and substitution warrants of its own
    /// `(read group, candidate)` parameter table, which `inference::repeat_tract_parameters`
    /// gathers — the calibration scale never enters a tract's likelihood, so charging a tract
    /// for it would report a worse warrant than the call has.
    ///
    /// **What is in neither is the prior's own shape**, which carries no provenance on the
    /// SNP/indel path at all; at a tract, the rung its shape came from travels on
    /// [`Self::repeat_tract`] rather than through this fold, because that is a
    /// statement about the prior where every rung here is a statement about the reads.
    pub weakest_provenance: Provenance,
    /// **What a repeat tract's call rested on, beyond the warrant above** — `None` at a SNP or
    /// indel, `Some` at every repeat tract.
    ///
    /// See [`RepeatTractProvenance`] for the four things it carries and why each of them cannot
    /// be recovered from the call.
    pub repeat_tract: Option<RepeatTractProvenance>,
    /// **How unlikely it is that no sample here carries a copy of any non-reference
    /// allele** — the site quality, on the Phred scale, as the worker computed it and
    /// *before* the artifact correction that the first output stage applies
    /// (`doc/devel/ng/spec/calling_quality.md` §5).
    ///
    /// **One quality field, written twice, and that is a rule with a shipped defect behind
    /// it.** Production keeps its corrected value nowhere and recomputes it at VCF-encode
    /// time; for sixteen days its `--min-qual` gate compared the engine's baseline while
    /// the *corrected* number went into the `QUAL` column, so sites were emitted `PASS`
    /// carrying a written quality of zero — 40 false positives at 30× on GIAB HG002 and 64
    /// at 50×, against 14 and 14 once both read one function (§3.5). ng carries one field:
    /// the worker writes the baseline here, the correction stage overwrites it in place,
    /// and there is never a second quality for anything to read by mistake.
    ///
    /// **Which is why this field is private and has no public reader.** §3.5's other half
    /// is that nothing between the worker and that stage may read it, and visibility is
    /// what enforces it: [`Self::uncorrected_site_quality`] is `pub(crate)`, so no consumer
    /// outside this crate can see the uncorrected number at all. The public accessor
    /// arrives with the stage that makes the value public-worthy. The check that the value
    /// came from the function that owns the ceiling is [`Self::new`]'s.
    site_quality: Phred,
    /// The nine pooled read counts the artifact correction consumes — or `None` at a locus
    /// that gives the two tests nothing to test.
    ///
    /// **Built here rather than downstream because one of the nine needs the calls**: how
    /// many alternative-allele reads the called genotypes lead you to expect. The other
    /// eight are the evidence's, and the evidence is released with the locus
    /// (`doc/devel/ng/spec/calling_quality.md` §3.3).
    ///
    /// **`None` is the ordinary outcome at more than one built locus in four**, not an
    /// error: the merge builds a locus when some sample's non-reference reads *pooled*
    /// reach its rule, so a locus whose candidate table came back as the reference alone is
    /// a first-class result — measured at 27.4% of built loci on the 63-accession tomato
    /// panel and 27.3% on HG002 at 30×
    /// ([`allele_candidates::SelectionVerdict::Selected`]). Both tests compare an
    /// alternative allele against the reference, so a locus with no alternative, or one
    /// whose alternatives drew no read at all, has nothing for them to weigh; production
    /// hands back its baseline unchanged in exactly those two cases
    /// ([`qual_refine.rs:79`](../../../src/vcf/qual_refine.rs)), and `None` is that same
    /// answer as a type rather than as an early return.
    artifact_test_counts: Option<ArtifactTestCounts>,
}

impl LocusInference {
    /// One locus's finished call, with everything checkable about it checked.
    ///
    /// # Panics
    ///
    /// If `cohort_expected_copies` is not one entry per allele of `alleles`, in either
    /// direction. [`ExpectedAlleleCopies`] is already built against *an* allele table, so
    /// what this catches is the residue: copies built against a **different** table that
    /// happened to be a different width. **Wider is the direction the pipeline
    /// produces** — the final prune shrinks the allele table, so copies not re-cut
    /// alongside it are left long, and every consumer indexes by [`AlleleId`], so the
    /// trailing entries would ride along unread. Both are in scope exactly here, which is
    /// why the check is here.
    ///
    /// If `per_sample` names no sample's call. A cohort has at least one sample and every
    /// sample gets a call, so an empty one is a locus that lost its calls rather than a
    /// locus with nothing to say.
    ///
    /// If `passes` is zero. Every locus takes at least one pass, so zero is a counter
    /// that was never incremented, and it would distort the distribution the pass cap is
    /// to be set from.
    ///
    /// If `repeat_tract` is set on a SNP/indel locus, **or missing at a repeat tract or
    /// bundle**. The
    /// tract ladder belongs to the repeat-tract prior; a generic locus seeds from a *frequency*
    /// spectrum, which has a ladder of its own with different rungs
    /// (`doc/devel/ng/spec/population_diversity.md` §3.4, §4.4), so a rung set there means one
    /// path's ladder was wired onto the other. The other direction is the newer half of the
    /// check: a tract whose record went missing reports no rung and no defaulted-cell counts,
    /// which reads downstream as *this call rested on nothing worth stating* rather than as a
    /// dropped field.
    ///
    /// If any call is [`SampleGenotypeCall::Missing`] at a repeat tract or bundle. **A
    /// repeat tract sets no sample aside** — a discovery round there can put back a length
    /// the cap cut, so no sample is locked out of the locus
    /// (`doc/devel/ng/spec/calling_em_loop.md` §5.0.1, §9). This is the mirror of the seed
    /// marker's check above: one ruling belongs to each path, and each is refused on the
    /// other. What it catches is a sample silently dropped from a tract's output — spec §9
    /// says such a sample has no expected copies at all, so a `Missing` wrongly emitted
    /// there also removes that sample from the cohort's expected-copies denominator, and
    /// the locus's allele frequencies come out wrong with a well-formed record beside them.
    ///
    /// If a called genotype names an allele this locus was not called over. `Genotype::new`
    /// checks only that the multiset is non-empty, and the final prune renumbers every id
    /// above the one it drops ([`CandidateAlleles::admit`]), so a call not remapped
    /// alongside it arrives stale. An out-of-range id reaches the VCF's `GT` column as an
    /// index past the `ALT` list — an unparseable record, the same failure class
    /// [`CandidateAlleles::new`] refuses an empty reference allele for. **It catches the
    /// out-of-range half only:** an id that stays in range after a renumber names a
    /// different allele silently, which is an argument for the prune returning its
    /// remapping rather than for a wider check here.
    ///
    /// If `region` runs backwards. [`GenomeRegion`] is plain data with no constructor of
    /// its own, and its documentation says a caller that requires `start <= end` says so
    /// itself; a called locus is such a caller, since the region reaches the `POS` column
    /// and the writer's span arithmetic.
    ///
    /// If `site_quality` is above [`quality::MAX_SITE_QUALITY`]. That ceiling belongs to
    /// [`quality::score_uncorrected_site_quality`], which is the only thing that may fill
    /// this field, so a value above it did not come from there — and the field is the one
    /// the correction stage later overwrites in place
    /// (`doc/devel/ng/spec/calling_quality.md` §3.5, §5.3). The lower end needs no check:
    /// [`Phred`] is finite and non-negative by construction.
    ///
    /// If an artifact summary names the reference as its primary alternative, or names an
    /// allele this locus was not called over. Both tests weigh one alternative against the
    /// reference, so a locus with nothing to weigh carries **no summary at all** rather
    /// than one pointing at allele 0; and an id past the table is the same stale-after-a-
    /// renumber failure the genotype check above catches, arriving by the other door.
    ///
    /// **What is deliberately *not* checked: `converged` against `passes`.** A locus that
    /// hit the cap should report the cap, but the cap is run configuration this type does
    /// not see. Nor is `converged` required to imply more than one pass: the loop's first
    /// pass compares against a reads-only estimate made before it
    /// (`doc/devel/ng/spec/calling_em_loop.md` §3), so settling on the first pass is a
    /// real outcome and not a comparison against nothing.
    #[allow(
        clippy::too_many_arguments,
        reason = "the architecture fixes this type's fields as a flat list \
                  (arch/calling_em_loop.md §2, and calling_quality.md §10's two additions); \
                  grouping them here to satisfy the lint would be a design change, not a \
                  refactor"
    )]
    pub fn new(
        region: GenomeRegion,
        alleles: CandidateAlleles,
        per_sample: Vec<SampleGenotypeCall>,
        cohort_expected_copies: ExpectedAlleleCopies,
        converged: bool,
        passes: u32,
        weakest_provenance: Provenance,
        repeat_tract: Option<RepeatTractProvenance>,
        site_quality: Phred,
        artifact_test_counts: Option<ArtifactTestCounts>,
    ) -> Self {
        assert_eq!(
            cohort_expected_copies.len(),
            alleles.len(),
            "a locus's expected allele copies run parallel to the alleles it was called \
             over: one entry per allele, reference first"
        );
        assert!(
            !per_sample.is_empty(),
            "a locus carries one call per sample and a cohort has at least one sample, so \
             a locus naming no call has lost them rather than having none to make"
        );
        assert!(
            passes > 0,
            "every locus takes at least one pass of the frequency loop, so a pass count \
             of zero is a counter that was never incremented"
        );
        // **A repeat *bundle* is on the demanding side of this, deliberately and unreachably.**
        // Written against `Generic` rather than against `Ssr(_)`, so a bundle would have to carry
        // a record — which it could not honestly build, having no single reference repeat count
        // and so no one allele length variation. Nothing can reach it: `tract_candidates` refuses a
        // bundle several frames earlier. When bundle scoring lands, this is one of the places
        // that has to answer for it, which is why it is stated rather than left to be discovered.
        assert_eq!(
            repeat_tract.is_some(),
            !matches!(alleles.kind(), LocusKind::Generic),
            "a {} locus {}: the tract ladder's rung and the defaulted-cell \
             counts belong to the repeat-tract path, and a SNP/indel locus is seeded from the \
             population curve's two moments, whose own ladder has different rungs",
            locus_kind_word(alleles.kind()),
            if repeat_tract.is_some() {
                "carries a repeat-tract record"
            } else {
                "carries no repeat-tract record"
            }
        );
        assert!(
            matches!(alleles.kind(), LocusKind::Generic)
                || per_sample.iter().all(|call| !call.is_missing()),
            "a repeat tract sets no sample aside, so a missing call at a {} locus is the \
             SNP/indel path's ruling wired onto the wrong path",
            locus_kind_word(alleles.kind())
        );
        assert!(
            per_sample
                .iter()
                .all(|call| call.genotype().is_none_or(|genotype| genotype
                    .alleles()
                    .iter()
                    .all(|id| alleles.bases_of(*id).is_some()))),
            "a called genotype names an allele this locus was not called over: the locus \
             holds {} alleles, and the final prune renumbers every id above the one it \
             drops, so a call not remapped alongside it is stale",
            alleles.len()
        );
        assert!(
            region.start <= region.end,
            "a called locus covers a stretch of reference, so its region cannot run \
             backwards: {region}"
        );
        assert!(
            site_quality.get() <= quality::MAX_SITE_QUALITY,
            "a site quality of {} is above the ceiling {} that \
             quality::score_uncorrected_site_quality caps at, so this number did not come \
             from it — and the field it is going into is the one the correction stage \
             overwrites in place",
            site_quality.get(),
            quality::MAX_SITE_QUALITY
        );
        if let Some(counts) = artifact_test_counts {
            assert!(
                !counts.primary_alternative.is_reference(),
                "the artifact tests compare an alternative allele against the reference, so \
                 a summary whose primary alternative *is* the reference has no comparison to \
                 make — a locus with nothing to test carries no summary at all rather than \
                 one naming allele 0"
            );
            assert!(
                alleles.bases_of(counts.primary_alternative).is_some(),
                "the artifact summary's primary alternative is allele {} and this locus \
                 holds {} alleles: the summary was pooled against a different allele table, \
                 or the final prune renumbered the table without re-cutting it",
                counts.primary_alternative.get(),
                alleles.len()
            );
        }
        Self {
            region,
            alleles,
            per_sample,
            cohort_expected_copies,
            converged,
            passes,
            weakest_provenance,
            repeat_tract,
            site_quality,
            artifact_test_counts,
        }
    }

    /// The alleles this locus was called over.
    #[inline]
    pub fn alleles(&self) -> &CandidateAlleles {
        &self.alleles
    }

    /// The cohort's expected allele copies, one entry per allele of [`Self::alleles`].
    #[inline]
    pub fn cohort_expected_copies(&self) -> &ExpectedAlleleCopies {
        &self.cohort_expected_copies
    }

    /// The site quality **as the worker wrote it**, before the artifact correction.
    ///
    /// **`pub(crate)`, and that is the enforcement of
    /// `doc/devel/ng/spec/calling_quality.md` §3.5 rather than a default choice of
    /// visibility**: the rule is that nothing between the worker and the correction stage
    /// may read this number, and the two that legitimately do — that stage, and this
    /// crate's tests — are both inside the crate. When the stage lands it overwrites the
    /// field in place and publishes a reader for the corrected value; there is no moment at
    /// which two qualities exist.
    ///
    /// **That stage now exists and this is read by it** (2026-09-01):
    /// [`evidence_for_output`](crate::run::records::evidence_for_output) takes this as the
    /// baseline, applies `correct_site_quality`, and the record carries the corrected number.
    /// The `dead_code` waiver this carried until then is spent.
    #[inline]
    #[must_use]
    pub(crate) fn uncorrected_site_quality(&self) -> Phred {
        self.site_quality
    }

    /// The nine pooled counts the artifact correction reads, or `None` where this locus
    /// gave its two tests nothing to weigh — see [`Self::artifact_test_counts`].
    #[inline]
    #[must_use]
    pub fn artifact_test_counts(&self) -> Option<ArtifactTestCounts> {
        self.artifact_test_counts
    }
}

#[cfg(test)]
mod tests;
