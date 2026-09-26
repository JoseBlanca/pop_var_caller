//! **The same cohort gives the same bytes on every machine**: a synthetic cohort is walked into
//! psps, its parameters are fitted, it is called with that fit, and the checksums of the
//! parameters file and of the VCF are compared against the ones recorded here.
//!
//! **Why the fit is in it.** Rust's `f64::ln`, `exp` and the rest call the operating system's
//! maths library, and glibc and Apple's library round differently in the last binary place. With
//! default parameters that almost never reaches a call, because a genotype is decided by a
//! comparison that one unit in the last place rarely flips: over 160 regions of four tomato
//! samples, macOS and Linux wrote identical VCFs even before the maths went through `libm`
//! (`doc/devel/reports/implementations/portable_float_B7_measurements_2026-09-15.md` §3.1). The
//! fit is different. It writes its numbers to the file in full, and the unchanged build's macOS
//! and Linux fits of the same psps differed in 7 of 574 lines. So a test that only called would
//! pass with the platform's library back in place; this one would not, and that was measured
//! (see `FITTED_PARAMETERS_MD5` below).
//!
//! **The fixture is the two-sample, 600-base varying cohort** (`cli::test_fixtures`): a SNP in
//! each sample, a repeat tract one sample shortens by a copy, two read groups in one sample. It
//! is one corner of the range this caller covers, chosen because the whole run takes under a
//! second; what it pins is that the arithmetic agrees, not that the calls are good.
//!
//! **What moves these checksums, and what to do when one moves.**
//!
//! - A change to the arithmetic of the walk, the fit or calling. That is what the test is for:
//!   the new checksums are recorded only after the change has been measured on the real cohort
//!   (`scripts/promote_ng_oracle.sh`) and its output difference explained.
//! - A different result on another platform, with the same source. That is the defect this test
//!   exists to catch. Do not re-record; find the call that rounds differently.
//! - **A new crate version**, which is written into the parameters file's census digests
//!   (`tool_version` in the repeat catalog's build settings) and into the VCF's `##source` line.
//!   Re-record both after checking that only those lines moved.
//!
//! **The pool's width is part of what is pinned.** A machine's core count decides how many
//! threads a run gets, so a second test runs the same cohort at one, four and seven threads and
//! requires the same bytes. It found a defect when it was written: the SNP/indel fit sized its
//! chunks from the pool's width and joined their floating-point totals in whatever order the pool
//! chose, and this fixture's fit came out as three different files at one, four and eight
//! threads (`tmp/digests_C/pool_width_cd056a3f.log`; the test itself uses seven). The fit's own
//! unit test `a_census_of_many_chunks_fits_to_the_same_bits_at_any_pool_width` checks the same
//! property on a census built to span several chunks, and was shown to fail with the old chunking.

#[cfg(test)]
mod tests {
    use std::path::{Path, PathBuf};

    use md5::{Digest, Md5};

    use crate::calling::allele_candidates::DEFAULT_MAX_CANDIDATE_ALLELES;
    use crate::cli::call_from_psps::{CallFromPspsArgs, run_call_from_psps};
    use crate::cli::estimate_parameters::{EstimateParametersArgs, run_estimate_parameters};
    use crate::cli::generate_psps::{GeneratePspsArgs, psp_path_for, run_generate_psps};
    use crate::cli::test_fixtures::{AVaryingCohort, a_varying_cohort_on_disk};
    use crate::region_typing::DEFAULT_MAX_STR_LEN;
    use crate::region_typing::segment_criteria::{
        DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
    };
    use crate::run::cohort_merge::DEFAULT_MAX_COHORT_LOCUS_SPAN;

    /// The width of the pool the pinned run uses; the other test checks other widths agree.
    const THREADS: usize = 4;

    /// The checksum of the parameters file the fit writes.
    ///
    /// Recorded 2026-09-15 with `float::exp` computed by table lookup
    /// (`doc/devel/implementation_plans/portable_float.md`, Milestone D), the same on macOS (arm64)
    /// and in the Linux container (arm64, glibc) (`tmp/digests_D1/{macos,linux}.log`); with
    /// `libm`'s `exp` it was `28722984e6600bd0cb5ba588506a9c81`. **Measured to have teeth:** the same test on the tree before the maths went through `libm`
    /// (`c6a4394b`, with the fit's fixed chunks and fixed join applied) wrote
    /// `f3f66a488983669e9e7e7a038338c010` on macOS and `cbcb75afaa4e1a25d5f314b581373028` on Linux
    /// (`tmp/digests_C/{macos,linux}_c6a4394b_final.log`).
    ///
    /// **Re-recorded 2026-09-25** when the SNP/indel fit gained SQUAREM acceleration and stopped at
    /// one part in a thousand rather than ten thousand — a deliberate change to the fit's
    /// arithmetic, measured on four tomato accessions before it was recorded (the fit ends 40
    /// log-likelihood units higher for the same passes). Recorded in the Linux container (arm64,
    /// glibc) only; it was `3edab375185d74ad84ab52255b418c0b`.
    const FITTED_PARAMETERS_MD5: &str = "70e634849f232de2fb32200dd21e0ab5";

    /// The checksum of the VCF called with that file, without its `##commandline` and
    /// `##reference` lines.
    ///
    /// Recorded with `libm`'s `exp`, and unchanged when `exp` moved to the table-driven version,
    /// which moved only the parameters file. **This one did not tell the platforms apart** before the
    /// maths went through `libm`: each platform's calls, made with its own differing fit, came to
    /// this same checksum. It is pinned because the calls are what the caller is for, and a change
    /// to them must be seen; the fit's checksum is the one guarding portability.
    ///
    /// **Re-recorded 2026-09-25** with the fit above, whose different numbers this file is called
    /// with; it was `3432b4219342c6277ec1d0072945e423`.
    const CALLS_MD5: &str = "dad5ff61fad91ca1d3d6a0e7587de946";

    /// Hex MD5 of some bytes.
    fn md5_hex(bytes: &[u8]) -> String {
        Md5::digest(bytes)
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect()
    }

    /// The VCF without the two header lines that name this run rather than its calls: the
    /// command line, which here is the test binary's own, and the reference's temporary path.
    fn comparable_vcf(vcf: &Path) -> String {
        std::fs::read_to_string(vcf)
            .expect("the VCF reads")
            .lines()
            .filter(|line| !line.starts_with("##commandline=") && !line.starts_with("##reference="))
            .map(|line| format!("{line}\n"))
            .collect()
    }

    /// Walk the cohort, fit it and call it with the fit. Hands back the parameters file's path
    /// and the VCF's.
    fn walked_fitted_and_called(cohort: &AVaryingCohort) -> (PathBuf, PathBuf) {
        let directory = cohort.directory.path();
        let psps = directory.join("psps");
        run_generate_psps(&GeneratePspsArgs {
            reference: cohort.reference.clone(),
            catalog: Some(cohort.catalog.clone()),
            alignments: cohort.alignments.clone(),
            output_dir: psps.clone(),
            regions: None,
            force: false,
            build_index_if_missing: false,
            min_copies: MinCopies::default(),
            min_period: DEFAULT_MIN_PERIOD,
            max_period: DEFAULT_MAX_PERIOD,
            max_str_len: DEFAULT_MAX_STR_LEN,
            min_purity: DEFAULT_MIN_PURITY,
        })
        .expect("the cohort walks into psps");

        let parameters = directory.join("fitted.parameters.toml");
        run_estimate_parameters(&EstimateParametersArgs {
            reference: cohort.reference.clone(),
            catalog: Some(cohort.catalog.clone()),
            psps: vec![psps.clone()],
            output: parameters.clone(),
            force: false,
            ploidy: 2,
            inbreeding: None,
            str_param_estimates_at_once: std::num::NonZeroUsize::MIN,
            skip_contamination: false,
        })
        .expect("the cohort fits");

        let calls = directory.join("calls.vcf");
        run_call_from_psps(&CallFromPspsArgs {
            reference: cohort.reference.clone(),
            catalog: Some(cohort.catalog.clone()),
            psps: ["one", "two"]
                .iter()
                .map(|sample| psp_path_for(&psps, sample))
                .collect(),
            output: calls.clone(),
            parameters: Some(parameters.clone()),
            defaults: false,
            ploidy: None,
            max_cohort_locus_span: DEFAULT_MAX_COHORT_LOCUS_SPAN,
            max_candidate_alleles: DEFAULT_MAX_CANDIDATE_ALLELES.get(),
            // The shipped false-discovery target, so the hidden-duplication filter scores and
            // writes its two fields: it computes with `ln`, `exp` and `ln_1p` as well.
            paralog_fdr: 0.01,
            paralog_filter_tag: false,
            cohort_locus_builder_regions_len: None,
            threads: 0,
            min_copies: MinCopies::default(),
            min_period: DEFAULT_MIN_PERIOD,
            max_period: DEFAULT_MAX_PERIOD,
            max_str_len: DEFAULT_MAX_STR_LEN,
            min_purity: DEFAULT_MIN_PURITY,
        })
        .expect("the cohort calls with its own fit");
        (parameters, calls)
    }

    /// One whole run — a fresh cohort walked, fitted and called — inside a pool of `threads`
    /// threads. Hands back the parameters file's text and the comparable VCF.
    fn one_run_at(threads: usize) -> (String, String) {
        let cohort = a_varying_cohort_on_disk();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("a pool of the asked-for width");
        let (parameters, calls) = pool.install(|| walked_fitted_and_called(&cohort));
        let fitted = std::fs::read_to_string(&parameters).expect("the fit wrote its file");
        (fitted, comparable_vcf(&calls))
    }

    /// **The fit and the calls made with it have the checksums recorded on the platforms this
    /// project runs on.**
    #[test]
    fn a_fitted_cohort_writes_the_same_bytes_on_every_platform() {
        let (fitted, called) = one_run_at(THREADS);
        let records = called.lines().filter(|line| !line.starts_with('#')).count();
        assert!(
            records >= 3,
            "the fixture varies at three positions, so a run that wrote {records} records is not \
             calling what this test means to pin"
        );
        assert!(
            called
                .lines()
                .filter(|line| !line.starts_with('#'))
                .any(|line| line.contains("PARALOG_POST=")),
            "the hidden-duplication filter must have scored the records"
        );

        let fitted_md5 = md5_hex(fitted.as_bytes());
        let calls_md5 = md5_hex(called.as_bytes());
        assert_eq!(
            (fitted_md5.as_str(), calls_md5.as_str()),
            (FITTED_PARAMETERS_MD5, CALLS_MD5),
            "the checksums moved; read this module's documentation before re-recording them.\n\
             --- fitted parameters ---\n{fitted}\n--- calls ---\n{called}",
        );
    }

    /// **The same bytes at one thread as at four and at seven**, so a machine's core count is not
    /// a reason for two fits of one cohort to differ.
    ///
    /// Before the SNP/indel fit's chunks were fixed in size and added in order, this fixture's fit
    /// wrote three different files at one, four and eight threads.
    #[test]
    fn a_fitted_cohort_writes_the_same_bytes_at_any_pool_width() {
        let at_four = one_run_at(THREADS);
        for threads in [1, 7] {
            let other = one_run_at(threads);
            assert!(
                other == at_four,
                "{threads} threads wrote a different fit or different calls from {THREADS}:\n\
                 --- fit at {threads} ---\n{}\n--- fit at {THREADS} ---\n{}",
                other.0,
                at_four.0,
            );
        }
    }
}
