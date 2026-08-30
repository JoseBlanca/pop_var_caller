//! **What the site quality costs across the cohort sizes this caller commits to.**
//!
//! `doc/devel/ng/spec/calling_quality.md` §13's Q3, and its own answer to how to settle it:
//! *"timing one locus's site quality at 63, 200, 1,000 and 3,000 samples — a benchmark over a
//! synthesised likelihood table, needing no data."* This is that benchmark.
//!
//! ```text
//! cargo bench --bench ng_site_quality_perf
//! ```
//!
//! # Why it is worth timing before anything calls a variant
//!
//! The site quality's second step folds the samples one at a time into a running distribution
//! over the **cohort's** total non-reference allele count, and that axis is
//! `ploidy × samples + 1` long. So the fold is quadratic in cohort size: 6,001 entries at 3,000
//! diploid samples, touched once per sample, which the spec prices at tens of millions of
//! multiply-adds a locus. **Production profiled this path at 200 samples and nowhere near the
//! top of the range** — the shipped caller's own comment records the linear-domain rewrite
//! winning back 88% of the path's time *at 200* — and nothing in this repository has run it
//! near the 3,000 `CLAUDE.md` §0 commits to.
//!
//! **If it is not affordable, the lever is the count axis and not the fold.** A cohort whose
//! fitted spectrum puts essentially no mass above a few dozen copies does not need 6,001
//! entries; truncating with a stated error bound is the honest version of that, and it is a
//! design change worth knowing about while the run is still being built rather than after it is
//! wired around the present shape. That is why this is timed now, months before ng emits a
//! file.
//!
//! # The two axes
//!
//! **`site_quality/samples`** — 63, 200, 1,000 and 3,000 diploid samples at a biallelic locus.
//! 63 is the tomato panel, 3,000 the top of the committed range. This is the axis the
//! quadratic lives on.
//!
//! **`site_quality/alleles`** — 2, 4 and 6 alleles at a fixed 1,000 samples. The *collapse*, the
//! only step that reads the whole `samples × genotypes` table, is linear in the genotype count,
//! which grows quadratically in the allele count: 3 genotypes at 2 alleles, 21 at 6. This axis
//! says whether the collapse or the fold is the half worth attacking, and the truncation lever
//! above only helps the second.
//!
//! # Two things in the timed body, both load-bearing
//!
//! **The scratch is prepared outside the timed closure.** `prepare_for_locus` sizes five buffers
//! and a run reuses one scratch across every locus of a segment, so timing the sizing would
//! measure an allocation the caller pays once and report it as per-locus work.
//!
//! **`black_box` on the returned quality.** The fold's result is otherwise dead in a benchmark,
//! and its four steps are pure arithmetic over buffers whose contents the optimiser can see.

use std::hint::black_box;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use pop_var_caller::ng::calling::genotype_prior::{SeedRegime, SpectrumSeed};
use pop_var_caller::ng::calling::quality::score_uncorrected_site_quality;
use pop_var_caller::ng::calling::{CallingScratch, CandidateAlleles, GenotypeTable};
use pop_var_caller::ng::locus_generation::LocusKind;
use pop_var_caller::ng::types::{LogProb, Ploidy};

/// The cohort sizes §13's Q3 names, plus nothing else: 63 is the tomato panel this repository
/// tests on, 3,000 the top of the range `CLAUDE.md` §0 commits to.
const COHORT_SIZES: [usize; 4] = [63, 200, 1_000, 3_000];

/// Allele counts for the second axis, at a fixed cohort. Six is the shipped candidate cap on
/// the SNP/indel path, so this spans what a run can actually present.
const ALLELE_COUNTS: [usize; 3] = [2, 4, 6];

/// The cohort the allele axis is swept at — mid-range, where the fold is already large enough
/// that a collapse cost has something to be compared against.
const ALLELE_AXIS_COHORT: usize = 1_000;

/// A neutral panel at one variant per kilobase — human diversity, and the middle column of
/// spec §5.4's table. The prior's two numbers do not touch the fold's cost; this is here so the
/// benchmark scores the arithmetic a run would score rather than a degenerate one.
fn human_like_seed() -> SpectrumSeed {
    SpectrumSeed::new(1.0, 1e-3, SeedRegime::NeutralShape)
}

fn diploid() -> Ploidy {
    Ploidy::try_new(2).expect("a diploid")
}

/// A candidate table of `alleles` one-base sequences, reference first.
fn candidates(alleles: usize) -> CandidateAlleles {
    const BASES: &[u8] = b"ACGTNM";
    let mut table = CandidateAlleles::new(Box::from(&BASES[0..1]), LocusKind::Generic);
    for allele in 1..alleles {
        table.admit(Box::from(&BASES[allele..allele + 1]));
    }
    assert_eq!(
        table.len(),
        alleles,
        "each fixture allele must be a distinct sequence, or the table collapses them"
    );
    table
}

/// **A cohort shaped like a rare variant**: one sample in fifty leans heterozygous, the rest
/// lean homozygous reference, and the leans differ from sample to sample.
///
/// The shape matters less than the fact that it is not uniform. A table of identical rows is
/// one the branch predictor and the rescaling both see coming, and the fold's per-sample
/// division by its own maximum is exactly the step a constant table would make free.
fn genotype_log_likelihoods(samples: usize, genotypes: usize) -> Vec<f64> {
    let mut table = Vec::with_capacity(samples * genotypes);
    for sample in 0..samples {
        let carrier = sample % 50 == 0;
        // A lean that walks, so no two samples hand the fold the same numbers.
        let lean = 4.0 + (sample % 7) as f64 * 0.6;
        for genotype in 0..genotypes {
            table.push(if carrier {
                // The first heterozygous genotype of the table is index 1 at every allele count.
                if genotype == 1 { 0.0 } else { -lean }
            } else if genotype == 0 {
                0.0
            } else {
                -lean * (1.0 + genotype as f64 / genotypes as f64)
            });
        }
    }
    table
}

/// One locus's scratch, prepared and filled — everything the timed body reads, built once.
fn prepared_scratch(
    samples: usize,
    alleles: usize,
) -> (
    CallingScratch<()>,
    CandidateAlleles,
    std::sync::Arc<GenotypeTable>,
) {
    let table = candidates(alleles);
    let genotypes = GenotypeTable::build(diploid(), alleles);
    let view = genotypes.view();
    let rows = genotype_log_likelihoods(samples, view.genotype_count());

    let mut scratch: CallingScratch<()> = CallingScratch::default();
    scratch.prepare_for_locus(samples, &table, &view);
    for sample in 0..samples {
        let row = &rows[sample * view.genotype_count()..(sample + 1) * view.genotype_count()];
        for (slot, &value) in scratch
            .sample_genotype_likelihoods_mut(sample)
            .iter_mut()
            .zip(row)
        {
            *slot = LogProb(value);
        }
    }
    let _ = view;
    (scratch, table, genotypes)
}

fn site_quality_by_cohort_size(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("site_quality/samples");
    for samples in COHORT_SIZES {
        let (mut scratch, _table, genotypes) = prepared_scratch(samples, 2);
        let view = genotypes.view();
        let seed = human_like_seed();
        group.bench_with_input(
            BenchmarkId::from_parameter(samples),
            &samples,
            |bencher, _| {
                bencher.iter(|| {
                    black_box(score_uncorrected_site_quality(
                        scratch.site_quality_buffers_mut(),
                        &view,
                        seed,
                    ))
                });
            },
        );
    }
    group.finish();
}

fn site_quality_by_allele_count(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("site_quality/alleles");
    for alleles in ALLELE_COUNTS {
        let (mut scratch, _table, genotypes) = prepared_scratch(ALLELE_AXIS_COHORT, alleles);
        let view = genotypes.view();
        let seed = human_like_seed();
        group.bench_with_input(
            BenchmarkId::from_parameter(alleles),
            &alleles,
            |bencher, _| {
                bencher.iter(|| {
                    black_box(score_uncorrected_site_quality(
                        scratch.site_quality_buffers_mut(),
                        &view,
                        seed,
                    ))
                });
            },
        );
    }
    group.finish();
}

criterion_group!(
    benches,
    site_quality_by_cohort_size,
    site_quality_by_allele_count
);
criterion_main!(benches);
