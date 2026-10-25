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
    use crate::cli::test_fixtures::{
        AVaryingCohort, a_larger_varying_cohort_with_sequencing_errors_on_disk,
        a_varying_cohort_on_disk, a_varying_cohort_with_sequencing_errors_on_disk,
    };
    use crate::region_typing::segment_criteria::{
        DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
    };
    use crate::region_typing::{DEFAULT_MAX_STR_BUNDLE_LEN, DEFAULT_MAX_STR_LEN};
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
    ///
    /// **Re-recorded 2026-09-28** when the SNP/indel fit's allele-frequency shapes began to step
    /// along the log-likelihood's own slope in them rather than solve the digamma form
    /// (`fit_precision.md` step A7), measured on four tomato accessions first
    /// (`scripts/promote_ng_oracle.sh`: the fit ends 30 log-likelihood units higher, at
    /// −7.991302 × 10⁵ against −7.991599 × 10⁵). A shape whose step would carry it past its bound is
    /// held there and the other takes the step the curvature gives with it held, so the free shape
    /// comes to rest where its own slope is zero. On this fixture two lines move: the ordinary-site
    /// prior's `reference_concentration`, 3.9591 to 39.089, and `alternative_concentration_total`,
    /// 7.30 × 10⁻¹⁰ to 6.60 × 10⁻³, both read off the fitted density, whose shapes moved from
    /// a = 10.8, b = 47.9 to a = 0.288, b = 50.0 (its upper bound); the fit's best log-likelihood
    /// rises from −44.919 to −44.704. Two samples over 600 bases do not place the shapes: the fit
    /// reports neither as identified. Recorded in the Linux container (arm64, glibc) only; it was
    /// `70e634849f232de2fb32200dd21e0ab5`. (Plan step A6, the read groups' own error rates, moved
    /// neither checksum: this fixture's reads carry no sequencing error, so every read group's two
    /// rates sit on their lower bounds either way.)
    ///
    /// **Re-recorded 2026-09-29** when the SNP/indel fit began to stop once every parameter is within
    /// a tenth of its standard error of the maximum, as a Newton step estimates the distance
    /// (`fit_precision.md` step B1), measured on four tomato accessions first
    /// (`scripts/promote_ng_oracle.sh`: every start there still runs to the 200-pass limit, and the
    /// fit's two checksums do not move). A parameter without an error is settled by definition, and
    /// two samples over 600 bases give almost none one — the mismapped share and the density's
    /// shapes are not identified, the invariant and fixed shares' errors are wider than their
    /// ranges — so every start now converges after 18 passes where it ran to the limit. Two lines
    /// move, both read off the density, whose first shape stops at 1.546 where 200 passes took it to
    /// 0.288 (the second at its bound, 50, both times): the ordinary-site prior's
    /// `reference_concentration`, 39.089 to 20.047, and `alternative_concentration_total`,
    /// 6.60 × 10⁻³ to 2.08 × 10⁻². The fit's best log-likelihood is −44.939 against −44.704.
    /// Recorded in the Linux container (arm64, glibc) only; it was
    /// `9d421df0bc755522f8a7b5617b4c30f8`.
    ///
    /// **Re-recorded 2026-10-06** when the parameters file became format version 2
    /// (`fit_precision.md` step E2): each number's standard error, how each start of the SNP/indel
    /// fit ended, and the notes explaining both. **No number moved**, measured first on four tomato
    /// accessions (`scripts/promote_ng_oracle.sh`: every call unchanged, both parameters files the
    /// old ones line for line once the version-2 keys and notes are taken out), and here: the same
    /// taking-out turns this fixture's new file back into the bytes of
    /// `22dda760d260572b786958e69b15b2b2`. The calls do not move. The file the review of the step
    /// saw (`f03c0ab14192957f7f1d8404bab4249a`) differs from this one in 19 lines of its notes and
    /// nothing else. Recorded in the Linux container (arm64, glibc) only; it was
    /// `22dda760d260572b786958e69b15b2b2`.
    ///
    /// **Re-recorded 2026-10-06** when a repeat-tract substitution rate stopped being zero (owner,
    /// checkpoint E of `fit_precision.md`): a count that found no mismatch takes half of one. This
    /// fixture's one stratum compared 460 bases and found none, so its three rows — one a read group
    /// — go from 0.0 to 0.5 / 461 = 0.0010846, now with a standard error of 0.0015347; nothing else
    /// moves but 8 lines of notes. Measured first on four tomato accessions: 52 of 188 rates move
    /// and nothing else in the file. Recorded in the Linux container (arm64, glibc), and the same on
    /// macOS (arm64); it was `a5e87936e968f9b855f91f88db828083`.
    ///
    /// **Re-recorded 2026-10-11** when the typing gained the STR bundle cap
    /// (`str_bundle_loci.md` §5.2). No number moves. Three things in the file do: `[repeat_routing]`
    /// gains a `max_str_bundle_len = 100` line, its note gains a clause on that key and now counts
    /// six flags of nine, and the census digests of the STR routing criteria and of the repeat
    /// catalog's build settings move, because each is taken over its criteria's printed form,
    /// which now carries the cap. Reverting those three gives the old checksum byte for byte; the
    /// calls do not move. Then again the same day, when the note's clause was reworded to say the
    /// hull runs to the furthest end of any repeat: restoring the earlier wording in the written
    /// file gives `4e643855f0dac8ee128f43a7a1a92095` exactly. Recorded in the Linux container (arm64, glibc) only; it was
    /// `860c3c65f27a1700e401d102833215e1`.
    ///
    /// **Re-recorded 2026-10-16** (`str_penalised_fit.md` plan step D4) after four changes, each
    /// measured on four tomato accessions first (`scripts/promote_ng_oracle.sh`;
    /// `doc/devel/reports/implementations/str_penalised_fit_D4_2026-10-16.md`):
    /// - "length spectrum" renamed "allele length variation" (`aa42d629`): two keys and four note
    ///   paragraphs. The allele window at ±9 (`38b1fe56`), merged at the same time, moved nothing
    ///   here: the climb fitted no row for this fixture's one tract.
    /// - file format 3 (plan steps B2, B4): `format_version = 3`, the fallback stated as the
    ///   between-locus share F (0.5, `defaulted`) in place of a concentration, and notes.
    /// - the repeat-tract fit is the penalised fit (C7), with a weak prior on stutter (C8): where the
    ///   climb wrote no repeat-tract row, the file now holds one slippage row and one allele-length
    ///   row for the two-base stratum at 10 repeats, with F 0.254, and the fallback F becomes that
    ///   F, `fitted_here`. Its slip rate is 0.0042: 5.8 × 10⁻⁷ at C7, before the prior's one read's
    ///   worth of default stutter.
    /// - note text (C8, C9, `7265efd0`).
    ///
    /// No number outside the repeat-tract rows moved; C10's ordering moved nothing, as the cohort's
    /// two samples were already in name order. Recorded on macOS (arm64) only; it was
    /// `a5b4997bc35e77425525e25e821cf40c`.
    ///
    /// **Re-recorded 2026-10-16 again** (plan step D8, where a thin class's errors start): two note
    /// paragraphs now state the lines of loci below which an error is withheld. No number or key
    /// moved: the contig holds one repeat tract ([`crate::cli::test_fixtures::TRACT`]), so the
    /// file's one repeat-tract row is a class of 1 locus with reads. That is below both lines of
    /// loci — 8 before D8, and 3 for the shorter share and fall-off after it — so every error of
    /// the row is withheld for too few loci, before its slipped reads are looked at, and the row
    /// carries no error before or after. Putting the two paragraphs' earlier wording back into the
    /// new file gives `b673dc9611c7126b5c35e9d7629493b3` exactly. The same holds for the next two
    /// fitted checksums, whose cohorts use the same contig. **So none of the three fitted
    /// checksums pins a repeat-tract error**, not even the two "with standard errors" pairs: their
    /// repeat-tract half pins values only. A cohort with at least 3 tracts of one class (8 for the
    /// slip rate) would be needed to put one in the file. Recorded on macOS (arm64) only.
    ///
    /// **Re-recorded 2026-10-16 a third time** (checkpoint E review, its nit on the notes): the
    /// repeat-tract section's notes explain a subset of samples (`samples_fitted_on`), numbers taken
    /// from a parametric curve or blended with one, and `curve_fitted_on` only when a row of the file
    /// holds one, which only a file written by the climb can. This file's one repeat-tract row is a
    /// point on a penalised curve, so its section note shrinks from 93 lines to 69. Note text
    /// alone: putting the old note block back into the new file gives the previous checksum
    /// exactly, for this and the next two fitted checksums
    /// (`tmp/review_E/fixes/digests/swap_notes.py`). Recorded on macOS (arm64) only; it was
    /// `168e7f326230e988286ef16b4c4fc5cd`.
    ///
    /// **Re-recorded 2026-10-16** (`snp_fit_convergence.md` plan step A2) when a start of the
    /// SNP/indel fit began to stop on the Newton decrement — the log-likelihood the remaining
    /// bounded Newton step is worth — below 0.001, in place of every parameter within a tenth of
    /// its error (`fit_precision.md` §2). Measured on four tomato accessions first
    /// (`scripts/promote_ng_oracle.sh`): every start there still runs to the limit, at the same pass
    /// counts, and the fitted file moves by its note on `snp_indel_fit_starts` alone. **Here every
    /// start runs to the limit again** (198 passes, where it converged at 18), and the file is back
    /// at the numbers of 2026-09-28: no parameter has an error on two samples over 600 bases, so the
    /// old test called every one settled at the first judging pass, while the decrement stays near
    /// 290 to the limit. That 290 is not what a step is worth: the bounded solve ends at its round
    /// limit with the density's second shape, already on its upper bound of 50, given a step of
    /// 2.7 × 10⁶, which the active-set rule returns only when it stops at that limit. Five lines
    /// move besides the note's: the three starts' rows (outcome and passes), and the ordinary-site
    /// prior's `reference_concentration`, 20.047 to 39.089, and
    /// `alternative_concentration_total`, 2.08 × 10⁻² to 6.60 × 10⁻³. Recorded on macOS (arm64)
    /// only; it was `66d72ee84c4a49d4b6271d3e5b68437c`.
    ///
    /// **Re-recorded 2026-10-16 again** (`snp_fit_convergence.md` plan step B1) when each judging
    /// pass that stops nothing began to try its Newton step, halved up to four times, each try a
    /// pass (`parameter_prepass_joint_fit.md` §3.3.1). Measured on four tomato accessions first: the
    /// built fit reproduces the prototype's run there and on one accession byte for byte but for the
    /// note. **Here every start still runs to the limit** (200 passes each, against 198): of the 29
    /// judging passes each start makes, 26 or 27 have every try refused — the quadratic model of
    /// the information on two samples over 600 bases is no guide — and each refusal spends up to
    /// four passes, so the starts reach the limit after fewer SQUAREM cycles and end elsewhere; 2 or 3
    /// steps were kept. Five lines move: the three starts' rows, and the ordinary-site prior's
    /// `reference_concentration`, 39.089 to 43.393, and `alternative_concentration_total`,
    /// 6.60 × 10⁻³ to 4.58 × 10⁻³. Recorded on macOS (arm64) only; it was
    /// `4966525dd98461b6aea6ecf845353b02`.
    ///
    /// **Re-recorded 2026-10-17** (`snp_fit_convergence.md` plan step B2) when the box solve behind
    /// the Newton step and its decrement became the primal active-set method (`fit_precision.md` §2,
    /// amendment of 2026-10-17). Measured on four tomato accessions first: the same passes, every
    /// SNP/indel number within 5 × 10⁻¹⁰ of itself, and one accession's file byte for byte. **Here
    /// every start now converges**, in 37, 41 and 41 passes, where all three ran to the limit: the
    /// old solve ended 164 of its 180 solves at its round limit and read the decrement off a step
    /// outside the box. The best log-likelihood rises from −44.7037 to −44.6991. Nine lines move:
    /// the starts' three rows; the three multipliers, 0.003164 to 0.003162 (every rate on its lower
    /// bound either way); both inbreeding coefficients, 3.2 × 10⁻¹⁰ to 0; the ordinary-site prior's
    /// `reference_concentration`, 43.393 to 48.274, and `alternative_concentration_total`,
    /// 4.58 × 10⁻³ to 1.79 × 10⁻²; and, in their tenth significant figure, the repeat-tract row and
    /// the fallback F that read the libraries' rates. Recorded on macOS (arm64) only; it was
    /// `5e36155136a9a0b33d5b89b6395d83fc`. (Plan step B3, the error-rate swap deleted, moved
    /// nothing here.)
    ///
    /// **Re-recorded 2026-10-17 again** (plan step B4) when the decrement began to count only the
    /// parameters that have a standard error, the others profiled out (`fit_precision.md` §2,
    /// amendment of 2026-10-17). Measured on four tomato accessions and one first: both files byte
    /// for byte, every parameter the step moves having an error there. **Here none of the 16 has
    /// one**, so the decrement is zero at the first judging pass and every start converges after 16
    /// passes (37, 41 and 41 before). The best log-likelihood is −44.871 against −44.699 at 41
    /// passes: the 0.17 units given up are all along parameters the data do not place. Nine lines
    /// move, back to the values of 2026-10-16 for six of them: the starts' rows; the multipliers
    /// 0.003162 to 0.003164 and the inbreeding coefficients 0 to 3.2 × 10⁻¹⁰; the ordinary-site
    /// prior's `reference_concentration`, 48.274 to 17.240, and `alternative_concentration_total`,
    /// 1.79 × 10⁻² to 1.63 × 10⁻²; the repeat-tract row and fallback F in their tenth significant
    /// figure. Recorded on macOS (arm64) only; it was `69445d42684b942826e43bdc02d49a3a`.
    ///
    /// **Re-recorded 2026-10-17 a third time** (plan step D2) for the parameters file's format 4:
    /// the version line, 3 to 4, and the two informational tables
    /// (`parameter_prepass_joint_fit.md` §2.1.3) — the density's, with its pair marked not
    /// separately determined (invariant share 0.971, first shape 1.33) and no number carrying an
    /// error, and the duplicated class's, marked the same way. With the two tables deleted and the
    /// version set back, the file's checksum is the old one; nothing a run reads moved, and the
    /// calls' checksum below did not. Recorded on macOS (arm64) only; it was
    /// `accb6b70b146bb0c83bfcedad8333871`.
    ///
    /// **Re-recorded 2026-10-17** when the tract ladder's stated concentration went from one to ten
    /// (`doc/devel/reports/ng_str_single_sample_prior_2026-10-17.md`): this fit has a fitted median,
    /// so the number is not read, and only the repeat-tract note naming the built-in value moved
    /// (0.5 to 1/11). Recorded on macOS (arm64) only; it was `66d72ee84c4a49d4b6271d3e5b68437c`.
    ///
    /// **Merged 2026-10-18 again**, main (the stated concentration above) into
    /// `snp-fit-convergence`: the notes from plan step A2 to D2 are the branch's, the one above
    /// main's. Against the branch's file only main's two note passages moved (the built-in 0.5
    /// became 1/11): putting the old two back gives the branch's `18ed96038479b4f9c860f1cb830038b0`
    /// exactly (`tmp/merge3/swap_back.py`). Recorded on the macOS host (arm64); it was
    /// `096c7151d55f21cd000c3209b07923e7` on main.
    ///
    /// **Re-recorded 2026-10-18** when the tract ladder's second rung got its own stated
    /// between-locus share (format 5, `doc/devel/reports/ng_str_single_sample_prior_2026-10-17.md`
    /// §8): this fixture's fit has a median of its own, so both keys carry it and no call moves.
    /// The file moves by `format_version`, the new key and the notes naming the two built-in
    /// values, and by nothing else. Recorded on the macOS host (arm64); it was `3ac3a2ae4c3e7e2fa38af1d6a3b827cc`.
    /// **Re-recorded 2026-10-18** when the built-in repeat-tract outlier weight went from 0.20 to
    /// 0.07 (parameters file format 6): the file moves by `format_version` and the weight's line,
    /// and by nothing else. Recorded on the macOS host (arm64); it was `ae20f59532f41838b4b63c0b2cf3fb5c`.
    /// **Re-recorded 2026-10-22** when the default base-quality multiplier moved from 1.0 to 2.0:
    /// no number moved, because every read group here is fitted; only the file's prose about a
    /// `defaulted` multiplier did. Putting the old prose back reproduces the old checksum,
    /// `467baf41a66b1053543571737bc7b48f`, exactly. Recorded in the Linux container (amd64).
    const FITTED_PARAMETERS_MD5: &str = "999e03bef3448832ba0285f0cc77013c";

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
    ///
    /// **Re-recorded 2026-09-28** with the fit above (step A7): the ordinary-site prior it carries
    /// moved, and with it three records' quality: the SNP at chrV:121 QUAL 480.2 to 541.9 and the
    /// other sample's GQ 65 to 73, the SNP at chrV:456 QUAL 481.0 to 542.7 and the other sample's GQ
    /// 65 to 73, the repeat tract at chrV:201 QUAL 0.0 to 7.4. No genotype moved. It was
    /// `dad5ff61fad91ca1d3d6a0e7587de946`.
    ///
    /// **Re-recorded 2026-09-29** with the fit above (step B1): the same three records' quality moves
    /// with the ordinary-site prior — the SNP at chrV:121 QUAL 541.9 to 549.5 and the other sample's
    /// GQ 73 to 71, the SNP at chrV:456 QUAL 542.7 to 550.3 and the other sample's GQ 73 to 71, the
    /// repeat tract at chrV:201 QUAL 7.4 to 14.3. No genotype moved. It was
    /// `7a570c442a4060fcd62f4af85b50a4fa`.
    ///
    /// **Re-recorded 2026-10-06** with the fit above, whose repeat-tract substitution rate is no
    /// longer zero: the repeat tract at chrV:201 moves, QUAL 14.3 to 14.9, AF 0.249885 to 0.249902,
    /// and one sample's GQ 33 to 34; its genotypes do not. The two SNPs do not move. Measured by
    /// putting the old rule back, which reproduces the old checksum exactly, and on four tomato
    /// accessions: 4 of 6,706 calls move, all repeat tracts, QUAL by at most 0.4, no genotype. The
    /// same in the Linux container and on macOS (both arm64). It was
    /// `2c0a89466d21c301a1704ed4e6cd0c84`.
    ///
    /// **Re-recorded 2026-10-07** when the strand and read-position penalty became four times a
    /// beta-binomial tail (`calling_quality.md` §6.2, §6.5), measured on GIAB at four depths and on
    /// 63 tomato accessions first. The fit does not move. The same three records and genotypes are
    /// written; the two SNPs' `SPPEN` rises and their `QUAL` falls by about 29 (chrV:121 549.5 to
    /// 520.7), because this fixture's reads lean one way at both alleles and the prior pulls the
    /// expected share a little toward even. It was `a2d17168b04f33a8f62a0a6cd2511df4`. Recorded in
    /// the Linux container (arm64, glibc) only.
    ///
    /// **Re-recorded 2026-10-10** when the hidden-duplication filter began scoring STR records from
    /// their read split as well as coverage (`hidden_paralog_filter.md` §3.2), measured on 63 tomato
    /// accessions first. The fit does not move. Only the STR record's `PARALOG_LR` moves (in the
    /// two-sample cohort, -1.11 to -19.24: its heterozygote's reads split 8 to 7, as a real
    /// heterozygote's do); no record, genotype or filter decision changes. It was `a44fc3523149780329f43fd951abdab7`.
    /// Recorded in the Linux container (arm64, glibc) only.
    ///
    /// **Re-recorded 2026-10-10, again,** when a repeat tract's site quality began to use the tract's own
    /// prior rather than the run's (`calling_quality.md` §8.1): the repeat tract at chrV:201 moves,
    /// QUAL 14.9 to 36.7; nothing else in the file does — the two SNPs keep the run's prior, the fit
    /// is unchanged, and no genotype moves. Measured by putting the old prior back at tracts, which
    /// reproduces the old checksum exactly. It was `2200dfaaa1148093cccf828021c5a281`. Recorded in
    /// the Linux container (arm64, glibc) only.
    ///
    /// **Re-recorded 2026-10-11** when STR bundles began to be called (`str_bundle_loci.md` §7.5):
    /// the header's `RU` and `PERIOD` became `Number=.`, `STR`'s description names bundles, and
    /// `REPCNB` is declared. Those four header lines are all that moved — putting the old four back
    /// reproduces the old checksum exactly; the fixture holds no bundle, and no record moves. It
    /// was `bc75ec4d5a3c67dcaa776f56746b2bdd`. Recorded in the Linux container (arm64, glibc) only.
    ///
    /// **Re-recorded 2026-10-16** with the fit above (plan step D4), which now holds a row for the
    /// repeat tract at chrV:201: that record's QUAL 36.7 to 90.1 (98.1 at C7, before the stutter
    /// prior), AF 0.249902 to 0.250006, the heterozygote's GQ 34 to 75 and the other sample's 45 to
    /// 47. No genotype moved; the two SNPs did not move. Recorded on macOS (arm64) only; it was
    /// `b1a1a907e4e7adf22fb62c56f3305258`.
    ///
    /// **Re-recorded 2026-10-16** with the fit above (`snp_fit_convergence.md` plan step A2): the
    /// two SNPs move with the ordinary-site prior, chrV:121 QUAL 520.7 to 513.1 and chrV:456 520.1
    /// to 512.5, the other sample's GQ 71 to 73 at both. No genotype moved, and the repeat tract
    /// did not move. Recorded on macOS (arm64) only; it was `dae7e77c280d7f240349cc8433207a81`.
    ///
    /// **Re-recorded 2026-10-16 again** with the fit above (plan step B1): the two SNPs' QUAL
    /// 513.1 to 511.1 and 512.5 to 510.5, the other sample's GQ 73 to 74 at both. No genotype moved;
    /// the repeat tract did not move. Recorded on macOS (arm64) only; it was
    /// `407693c2f919227045b95996191f12ed`.
    ///
    /// **Re-recorded 2026-10-17** with the fit above (plan step B2): the two SNPs' QUAL, chrV:121
    /// 511.1 to 516.6 and chrV:456 510.5 to 516.0. No genotype or GQ moved; the repeat tract did
    /// not move. Recorded on macOS (arm64) only; it was `24afa05244c632f84f33f05bf86b4675`.
    ///
    /// **Re-recorded 2026-10-17 again** with the fit above (plan step B4): the two SNPs' QUAL
    /// 516.6 to 520.2 and 516.0 to 519.6, the other sample's GQ 74 to 70 at both. No genotype
    /// moved; the repeat tract did not move. Recorded on macOS (arm64) only; it was
    /// `8babc6961203813baa8802274c79e728`.
    ///
    /// **Re-recorded 2026-10-16** when the hidden-duplication filter's VCF id became `hiddenDup`
    /// (`doc/devel/ng/spec/command_line.md` §4.3): its `##FILTER` header line is all that moved —
    /// putting `hiddenParalog` back reproduces the old checksum exactly. It was `dae7e77c280d7f240349cc8433207a81`, the value `main` recorded for the fit above.
    /// Recorded on the macOS host (arm64).
    ///
    /// **Re-recorded 2026-10-16** when a calling run stopped writing a parameters file beside its
    /// VCF (`doc/devel/ng/spec/command_line.md` §7): `##source` reads `pvc` where it read `ng`, and
    /// `##parametersFile=calls.parameters.toml` became `##pvcLog` and `##pvcParameters`. Those
    /// header lines are all that moved — putting the old two back reproduces the old checksum
    /// exactly. It was `eb43b6968066b5d22525aaaeaed2f736`. Recorded on the macOS host (arm64).
    ///
    /// **Re-recorded 2026-10-17** when the coverage histogram's depth axis came to be fitted from
    /// every window through a fixed relative-depth grid (`doc/devel/ng/spec/window_coverage.md`
    /// §3.4): `PARALOG_LR` and `PARALOG_POST` are all that moved — with both taken out, the old
    /// and new files are identical line for line, the filter's header line included, and all
    /// 3 records are kept by both. The fixture is under 10,000 windows, so the old axis was
    /// already fitted from every window, and the move is the grid's rounding of the median and the
    /// re-cut alone. It was `e9738457fb6b4e36f33e5fa008b935f0`. Recorded on the macOS host (arm64).
    ///
    /// **Merged 2026-10-18**, main into `snp-fit-convergence`: the notes from plan step A2 to here
    /// are the branch's and those from the `hiddenDup` one on are main's, each side moving the
    /// checksum for its own reasons. Against the branch's VCF the merged one differs only by main's
    /// three changes — the `##source`/`##pvcLog`/`##pvcParameters` header lines, the filter's id
    /// `hiddenDup`, and `PARALOG_LR`/`PARALOG_POST` from the coverage grid; with those taken out the
    /// two files are identical, all 3 records included. The fitted file did not move. It was
    /// `2661c66dac4f40502bda5f7e2e242713` on the branch and `ce76015b7bc2f654f0adf2c235062381` on
    /// main. Recorded on the macOS host (arm64).
    ///
    /// **Re-recorded 2026-10-18**, merging the tract ladder's stated concentration (one to ten,
    /// `doc/devel/reports/ng_str_single_sample_prior_2026-10-17.md`) onto the above: only the
    /// `##pvcParameters` md5 moved — putting the old one back reproduces `ce76015b7bc2f654f0adf2c235062381` exactly.
    /// Recorded on the macOS host (arm64).
    ///
    /// **Merged 2026-10-18 again**, as [`FITTED_PARAMETERS_MD5`]'s last note: only the
    /// `##pvcParameters` md5 moved, with the fitted file's note — putting the branch's md5 back
    /// gives the branch's `7672e055a1dbe6a9a642fc8badba944c` exactly; no record moved. It was
    /// `bd3ec1e0648501deaff70599412f80a4` on main. Recorded on the macOS host (arm64).
    ///
    /// **Re-recorded 2026-10-18** with the fit above: only the `##pvcParameters` md5 moved —
    /// putting the old one back reproduces the old checksum exactly. It was `2b5fedb559b16c8dfed07ab39f7cc5d2`. Recorded on
    /// the macOS host (arm64).
    /// **Re-recorded 2026-10-18** with the outlier weight above: besides the `##pvcParameters` md5,
    /// the one STR record (`chrV:201`) moves in QUAL and in its heterozygous sample's GQ — more
    /// confident at the lower floor — with every genotype unchanged. Recorded on the macOS host
    /// (arm64); it was `d90fefc996dca55d37c45c43cd392170`.
    /// **Re-recorded 2026-10-21** when the hidden-duplication filter's allele term became a
    /// beta-binomial (`DEFAULT_ALLELE_DISPERSION`, `doc/devel/reports/ng_dup_filter_allele_balance_2026-10-21.md`):
    /// only the two records' `PARALOG_LR` (chrV:121 −21.66 → −21.62, chrV:201 −19.13 → −18.19) moved; every record still passes, and `PARALOG_POST`, the header and every
    /// other column are unchanged. Setting the dispersion to zero reproduces `0999e9cf0affd975c1ef46d6d829f923` exactly.
    /// Recorded on the macOS host (arm64).
    /// **Re-recorded 2026-10-22** when the repeat-tract periodicity check was removed
    /// (`doc/devel/ng/spec/candidate_alleles_ssr.md` §7): only the header's
    /// `##FILTER=<ID=notPeriodic,…>` line went; putting it back reproduces the old checksum,
    /// `cea47374a0885f91c92940f2d89e103e`, exactly. Recorded in the Linux container (amd64).
    /// **Re-recorded 2026-10-22** when the default base-quality multiplier moved from 1.0 to 2.0:
    /// only the `##pvcParameters` line moved, carrying the parameters file's new checksum; with
    /// the old one put back it reproduces the old checksum, `da2024ac5ff5582c96a0f86ce9a8254d`, exactly. Recorded in the
    /// Linux container (amd64).
    /// **Re-recorded 2026-10-23** when `PARALOG_POST`'s header description became "under the paralog
    /// rate fitted around it" for the windowed duplication rate: only that `##INFO` line moved; putting
    /// the old description back reproduces the old checksum, `a98b088ade1653000a043072240ad900`, exactly. Recorded
    /// in the Linux container (amd64).
    const CALLS_MD5: &str = "97c0120bca8f944d9d26e07695bf87c2";

    /// **The checksum of the parameters file fitted on the same cohort with sequencing errors in its
    /// reads** ([`a_varying_cohort_with_sequencing_errors_on_disk`]) — the file that carries the fit's
    /// standard errors, which the file above does not: there every library's rate sits on its
    /// bound. Added at checkpoint E of `fit_precision.md` (owner, 2026-10-06), so that the errors'
    /// arithmetic — the information summed over positions, inverted, carried to a multiplier — is
    /// pinned across platforms and pool widths as the numbers themselves are.
    ///
    /// **What it carries**: each library's multiplier, 2.97 to 6.43 (an error in about one read in
    /// four over reads claiming Q35), with a standard error of 1.42 to 17.71 — wide, since 600
    /// positions over two samples barely place a library's rate, but computed, which is what is
    /// pinned; and the tract's substitution rate with its binomial error. The two inbreeding
    /// coefficients carry none: one variant a sample does not place them. The calls are the
    /// cohort's three designed records. **Its errors come from the whole matrix**, as every cohort
    /// of at most 20 samples' do; the next pair pins the other path. **Not every row of that matrix
    /// reaches the file**: measured in this step's review, the cohort's score for the share of
    /// mismapped positions can be computed differently, even doubled, without moving a byte here, so
    /// a change confined to the cohort's own scores is not pinned by this pair. Recorded in the Linux container
    /// (arm64, glibc) and the same on macOS (arm64).
    ///
    /// **Re-recorded 2026-10-11, twice,** for the two reasons [`FITTED_PARAMETERS_MD5`] gives.
    /// Recorded in the Linux container (arm64, glibc) only; it was `8c71a5a95f1196e3dd8c05abdccf7355`, then `9265d0273c5400d8a2b2554b3f7ecf0a`.
    ///
    /// **Re-recorded 2026-10-16** for the four reasons [`FITTED_PARAMETERS_MD5`] gives, with the
    /// same repeat-tract rows to the last digit. Recorded on macOS (arm64) only; it was
    /// `4859eacd5645e607b6f99b667d184dc0`. **Again the same day** for the note text of plan step
    /// D8 alone, as [`FITTED_PARAMETERS_MD5`] says; it was `3a016fb7a40bc91eb5bb9ae2adc6e286`.
    /// **A third time** for the conditional notes alone, as [`FITTED_PARAMETERS_MD5`] says; it was
    /// `a48841ff7b0d8ca924b758dc7a33fcbf`.
    ///
    /// **Re-recorded 2026-10-16** (`snp_fit_convergence.md` plan step A2), for the reason
    /// [`FITTED_PARAMETERS_MD5`] gives: the three starts, which converged at 49, 54 and 15 passes,
    /// now run to the limit (199, 198, 198). Start 1 had been called converged at its judging pass 47,
    /// every parameter with an error within a tenth of it, while the decrement there was 290 — again
    /// a bounded solve ended at its round limit; it then climbs from −489.0 to −482.4 by pass 197.
    /// Its judging passes alternate between about −482.4 and −532.0 from pass 53 to the limit:
    /// every other cycle loses some 50 units, which the old rule's early stop had hidden; where in
    /// the cycle they are lost was not traced. The multipliers move
    /// from 2.97, 6.43 and 3.13 to 0.84, 6.37 and 1.09 (errors 0.54, 16.7 and 0.53), and the
    /// ordinary-site prior's two concentrations with the density; plus the note's seven lines.
    /// Recorded on macOS (arm64) only; it was `0586cd65f3dcbd192f0f8cf19e223ca9`.
    ///
    /// **Re-recorded 2026-10-16 again** (plan step B1), for the reason [`FITTED_PARAMETERS_MD5`]
    /// gives: every start still at the limit (198, 199, 200 passes), with 7, 16 and 3 Newton steps
    /// kept and 21, 23 and 27 judging passes whose every try was refused. The multipliers move from
    /// 0.84, 6.37 and 1.09 to 6.68, 6.24 and 6.09 (errors 13.4, 42.6 and 3.8), the ordinary-site
    /// prior's two concentrations with the density, and the starts' three rows. Recorded on macOS
    /// (arm64) only; it was `4c64cac0190ef850c66dad4b50d366f1`.
    ///
    /// **Re-recorded 2026-10-17** (plan step B2), for the reason [`FITTED_PARAMETERS_MD5`] gives.
    /// Here every start still runs to the limit (199, 199, 198 passes): with the solve finishing, the
    /// decrement is honest but stays above the threshold along parameters that have no error, and
    /// the error-rate swap still makes plain passes lose tens of log-likelihood units (both
    /// fixed in the next steps). The best log-likelihood rises from −483.91 to −482.16. The
    /// multipliers move from 6.68, 6.24 and 6.09 to 7.75, 6.44 and 8.38 (errors 7.4, 28.8 and 4.2),
    /// the inbreeding coefficients from 3.2 × 10⁻¹⁰ to 3.1 × 10⁻¹³, the ordinary-site prior's two
    /// concentrations with the density (34.85 to 49.89, 0.0191 to 0.0223), the starts' rows, and in
    /// the tenth significant figure the repeat-tract row and fallback F. Recorded on macOS (arm64)
    /// only; it was `e5463fbf4a49798b810cc44db00f3a6d`.
    ///
    /// **Re-recorded 2026-10-17 again** (plan step B3) when the maximisation stopped swapping a
    /// library's two error rates whenever the mismapped one came out below the ordinary one
    /// (`parameter_prepass_joint_fit.md` §3.3, decided 2026-10-17). Measured on four tomato
    /// accessions and one first: both files byte for byte, the swap never firing there. Here it
    /// fired on every other cycle and lost about 50 units each time; without it the three starts end
    /// within 0.001 of each other (−482.158 to −482.159) where they ended 1.0 unit apart, and two of
    /// the three libraries end with the mismapped rate below the ordinary one, which the run log now
    /// says. Every start still runs to the limit (198, 199, 199): the decrement stays above the
    /// threshold along parameters that have no error, which the next step addresses. The multipliers
    /// move by at most 0.4% (7.754, 6.411, 8.379), the inbreeding coefficients 3.1 × 10⁻¹³ to
    /// 1.1 × 10⁻¹¹, the ordinary-site prior's concentrations 49.89 to 47.53 and 0.0223 to 0.0194,
    /// the starts' rows, and the repeat-tract row and fallback F in their eleventh significant
    /// figure. Recorded on macOS (arm64) only; it was `64d573193f0492bd973d8d0a82367448`.
    ///
    /// **Re-recorded 2026-10-17 a third time** (plan step B4), for the reason
    /// [`FITTED_PARAMETERS_MD5`] gives. Here 12 of the 16 parameters have no error; measured over
    /// the other four the decrement falls below the threshold, and **every start converges**, in
    /// 22, 13 and 13 passes where all three ran to the limit. The best log-likelihood is −482.266
    /// against −482.158 after 198 passes. The multipliers move by at most 0.5% (7.753, 6.444,
    /// 8.378; errors 6.9, 27.3 and 3.9), the inbreeding coefficients 1.1 × 10⁻¹¹ to 3.2 × 10⁻¹⁰,
    /// the ordinary-site prior's concentrations 47.53 to 18.94 and 0.0194 to 0.0148, the starts'
    /// rows, and the repeat-tract row and fallback F in their tenth significant figure. Recorded
    /// on macOS (arm64) only; it was `d5a9d68eb60e8ac53f61fa1601fb7d50`.
    ///
    /// **Re-recorded 2026-10-17 a fourth time** (plan step D2), for the reason
    /// [`FITTED_PARAMETERS_MD5`] gives: the version line and the two informational tables, both
    /// pairs marked (invariant share 0.974, first shape 1.37; duplicated share 0.0041, carrier
    /// shape 3.10), no number with an error. The tables deleted and the version set back give the
    /// old checksum; the calls' did not move. Recorded on macOS (arm64) only; it was
    /// `16cdb5d1f7ab2ece9d505bdba3586eb4`.
    ///
    /// **Re-recorded 2026-10-17** when the tract ladder's stated concentration went from one to ten
    /// (`doc/devel/reports/ng_str_single_sample_prior_2026-10-17.md`): this fit has a fitted median,
    /// so the number is not read, and only the repeat-tract note naming the built-in value moved
    /// (0.5 to 1/11). Recorded on macOS (arm64) only; it was `0586cd65f3dcbd192f0f8cf19e223ca9`.
    ///
    /// **Merged 2026-10-18 again**, as [`FITTED_PARAMETERS_MD5`]'s last note: main's two note
    /// passages alone; putting the old two back gives the branch's
    /// `0138ccccf14e750e671290de06b58afd` exactly. Recorded on the macOS host (arm64); it was
    /// `af39cc216bf85d18914a6112e8403d5b` on main.
    ///
    /// **Re-recorded 2026-10-18** when the tract ladder's second rung got its own stated
    /// between-locus share (format 5, `doc/devel/reports/ng_str_single_sample_prior_2026-10-17.md`
    /// §8): this fixture's fit has a median of its own, so both keys carry it and no call moves.
    /// The file moves by `format_version`, the new key and the notes naming the two built-in
    /// values, and by nothing else. Recorded on the macOS host (arm64); it was `49bd109f04196c5d2a3f39341002f691`.
    /// **Re-recorded 2026-10-18** when the built-in repeat-tract outlier weight went from 0.20 to
    /// 0.07 (parameters file format 6): the file moves by `format_version` and the weight's line,
    /// and by nothing else. Recorded on the macOS host (arm64); it was `ca1be8e12b7f816755ba43d1e5fcd4c9`.
    /// **Re-recorded 2026-10-22** when the default base-quality multiplier moved from 1.0 to 2.0:
    /// no number moved, because every read group here is fitted; only the file's prose about a
    /// `defaulted` multiplier did. Putting the old prose back reproduces the old checksum,
    /// `d80331a16f4cd86bfbdf885c90cc7188`, exactly. Recorded in the Linux container (amd64).
    const WITH_ERRORS_FITTED_PARAMETERS_MD5: &str = "74601f6755190f65c159a41c95bff1ed";

    /// The checksum of the VCF called with that file, as [`CALLS_MD5`] is of the first.
    ///
    /// **Re-recorded 2026-10-07** with the strand penalty of [`CALLS_MD5`]'s last note: the same
    /// records and genotypes, the SNPs' `QUAL` about 29 lower. It was
    /// `943c7fa93e2a1f110e996c117e3b0808`. Recorded in the Linux container (arm64, glibc) only.
    ///
    /// **Re-recorded 2026-10-10** when the hidden-duplication filter began scoring STR records from
    /// their read split as well as coverage (`hidden_paralog_filter.md` §3.2), measured on 63 tomato
    /// accessions first. The fit does not move. Only the STR record's `PARALOG_LR` moves (in the
    /// two-sample cohort, -1.11 to -19.24: its heterozygote's reads split 8 to 7, as a real
    /// heterozygote's do); no record, genotype or filter decision changes. It was `b82fed4ad01da17688bc4d695a9ed03d`.
    /// Recorded in the Linux container (arm64, glibc) only.
    ///
    /// **Re-recorded 2026-10-10, again,** with the tract prior of [`CALLS_MD5`]'s last note: only the
    /// repeat tract's `QUAL` moves. It was `206e0f7bf9b5459b99b20477562642a1`. Recorded in the
    /// Linux container (arm64, glibc) only.
    ///
    /// **Re-recorded 2026-10-11** when STR bundles began to be called (`str_bundle_loci.md` §7.5):
    /// the header's `RU` and `PERIOD` became `Number=.`, `STR`'s description names bundles, and
    /// `REPCNB` is declared. Those four header lines are all that moved — putting the old four back
    /// reproduces the old checksum exactly; the fixture holds no bundle, and no record moves. It
    /// was `f037a347e52634175a3eb5b4353fb8d0`. Recorded in the Linux container (arm64, glibc) only.
    ///
    /// **Re-recorded 2026-10-12** when a partial read that held the right flank of an STR segment
    /// began to be compared with the candidate's end rather than its start (`read_likelihoods.md`
    /// §5.2). Only the repeat tract's `AF` moves, in its sixth decimal (0.249923 to 0.249924).
    /// Measured by scoring every partial from the left again, which restores the old checksum; the
    /// benchmark measurement is `reports/implementations/str_partial_read_phase_2026-10-12.md`.
    /// It was `f037a347e52634175a3eb5b4353fb8d0`. Recorded in the Linux container (arm64, glibc)
    /// only.
    ///
    /// **Merged 2026-10-12** (main's partial-read phase fix into the STR-bundle branch): the two
    /// notes above were each recorded from `f037a347e52634175a3eb5b4353fb8d0` on their own branch.
    /// With both, the header is the bundle branch's and the repeat tract's record is main's to the
    /// byte (`AF=0.249924`); segments now also record a whole-tract partial's held flank, and no
    /// field of this fixture moves for it. It was `977f08d389bdba037ee350ab38164e1b` on str-bundles
    /// and `e9e8aee5b6714aa33ce91bd04290984e` on main. Recorded in the Linux container (arm64,
    /// glibc) only.
    ///
    /// **Re-recorded 2026-10-16** with the fit above (plan step D4): only the repeat tract at
    /// chrV:201 moves, QUAL 37.4 to 90.4, AF 0.249924 to 0.250010, the heterozygote's GQ 34 to 74
    /// and the other sample's 42 to 44; no genotype. Recorded on macOS (arm64) only; it was
    /// `f65f586f58163d23d127d7f8c823e433`.
    ///
    /// **Re-recorded 2026-10-16** with the fit above (`snp_fit_convergence.md` plan step A2): the
    /// two SNPs' QUAL, chrV:121 222.1 to 277.1 and chrV:456 219.2 to 265.1, and the other sample's GQ
    /// 67 to 71 at both. No genotype moved; the repeat tract did not move. Recorded on macOS (arm64)
    /// only; it was `d315b4e4d0de2f21fd34cd673bbdfa8b`.
    ///
    /// **Re-recorded 2026-10-16 again** with the fit above (plan step B1): the two SNPs' QUAL
    /// 277.1 to 185.7 and 265.1 to 189.1, with the multipliers' move, and the other sample's GQ 71
    /// to 73 at both. No genotype moved; the repeat tract did not move. Recorded on macOS (arm64)
    /// only; it was `3ccf7a635b25ef0dcbfc54bc2fb42e92`.
    ///
    /// **Re-recorded 2026-10-17** with the fit above (plan step B2): the two SNPs' QUAL 185.7 to
    /// 178.5 and 189.1 to 174.5, with the multipliers' move, and the other sample's GQ 73 to 74 at
    /// both. No genotype moved; the repeat tract did not move. Recorded on macOS (arm64) only; it
    /// was `46d43b7dc786bc8717807898ebd67ac5`.
    ///
    /// **Re-recorded 2026-10-17 again** with the fit above (plan step B3): the two SNPs' QUAL
    /// 178.5 to 178.1 and 174.5 to 174.1. Nothing else moved. Recorded on macOS (arm64) only; it was
    /// `37d84c3ef1b0fe88160b08eeb72f1f90`.
    ///
    /// **Re-recorded 2026-10-17 a third time** with the fit above (plan step B4): the two SNPs'
    /// QUAL 178.1 to 180.5 and 174.1 to 176.5, the other sample's GQ 74 to 70 at both. No genotype
    /// moved; the repeat tract did not move. Recorded on macOS (arm64) only; it was
    /// `241e57ecda2afed9a792607a44a8bb8c`.
    ///
    /// **Re-recorded 2026-10-16** when the hidden-duplication filter's VCF id became `hiddenDup`
    /// (`doc/devel/ng/spec/command_line.md` §4.3): its `##FILTER` header line is all that moved —
    /// putting `hiddenParalog` back reproduces the old checksum exactly. It was `d315b4e4d0de2f21fd34cd673bbdfa8b`, the value `main` recorded for the fit above.
    /// Recorded on the macOS host (arm64).
    ///
    /// **Re-recorded 2026-10-16** when a calling run stopped writing a parameters file beside its
    /// VCF (`doc/devel/ng/spec/command_line.md` §7): `##source` reads `pvc` where it read `ng`, and
    /// `##parametersFile=calls.parameters.toml` became `##pvcLog` and `##pvcParameters`. Those
    /// header lines are all that moved — putting the old two back reproduces the old checksum
    /// exactly. It was `3a2e830acc6991ee3b2285a80abbb89e`. Recorded on the macOS host (arm64).
    ///
    /// **Re-recorded 2026-10-17** when the coverage histogram's depth axis came to be fitted from
    /// every window through a fixed relative-depth grid (`doc/devel/ng/spec/window_coverage.md`
    /// §3.4): `PARALOG_LR` and `PARALOG_POST` are all that moved — with both taken out, the old
    /// and new files are identical line for line, the filter's header line included, and all
    /// 3 records are kept by both. The fixture is under 10,000 windows, so the old axis was
    /// already fitted from every window, and the move is the grid's rounding of the median and the
    /// re-cut alone. It was `973688dbec909bb94ee73cfd2941401c`. Recorded on the macOS host (arm64).
    ///
    /// **Merged 2026-10-18**, main into `snp-fit-convergence`, as [`CALLS_MD5`]'s last note: only
    /// main's header lines, filter id and `PARALOG_LR`/`PARALOG_POST` differ from the branch's VCF.
    /// It was `25a1e5728030f1b0dbdfa916fbf4ee74` on the branch and
    /// `b1ae11de0ce915350cd2b76b2b72d630` on main. Recorded on the macOS host (arm64).
    ///
    /// **Re-recorded 2026-10-18**, merging the tract ladder's stated concentration (one to ten,
    /// `doc/devel/reports/ng_str_single_sample_prior_2026-10-17.md`) onto the above: only the
    /// `##pvcParameters` md5 moved — putting the old one back reproduces `b1ae11de0ce915350cd2b76b2b72d630` exactly.
    /// Recorded on the macOS host (arm64).
    ///
    /// **Merged 2026-10-18 again**, as [`CALLS_MD5`]'s last note: the `##pvcParameters` md5 alone;
    /// putting the branch's back gives the branch's `c81df78823bf9aa0836eee9dca56e85f` exactly. It
    /// was `69eebab6d26ad5a44d9512631851fbba` on main. Recorded on the macOS host (arm64).
    ///
    /// **Re-recorded 2026-10-18** with the fit above: only the `##pvcParameters` md5 moved —
    /// putting the old one back reproduces the old checksum exactly. It was `23c98fcb501498aa7167d74fedcd51fe`. Recorded on
    /// the macOS host (arm64).
    /// **Re-recorded 2026-10-18** with the outlier weight above: besides the `##pvcParameters` md5,
    /// the one STR record (`chrV:201`) moves in QUAL and in its heterozygous sample's GQ — more
    /// confident at the lower floor — with every genotype unchanged. Recorded on the macOS host
    /// (arm64); it was `38ee2760d8b096add2b6019fa58765b5`.
    /// **Re-recorded 2026-10-21** when the hidden-duplication filter's allele term became a
    /// beta-binomial (`DEFAULT_ALLELE_DISPERSION`, `doc/devel/reports/ng_dup_filter_allele_balance_2026-10-21.md`):
    /// only the two records' `PARALOG_LR` (chrV:121 −21.66 → −21.62, chrV:201 −19.13 → −18.19) moved; every record still passes, and `PARALOG_POST`, the header and every
    /// other column are unchanged. Setting the dispersion to zero reproduces `a7215a1ca6142fa3ebe3e846ef96bd24` exactly.
    /// Recorded on the macOS host (arm64).
    /// **Re-recorded 2026-10-22** when the repeat-tract periodicity check was removed
    /// (`doc/devel/ng/spec/candidate_alleles_ssr.md` §7): only the header's
    /// `##FILTER=<ID=notPeriodic,…>` line went; putting it back reproduces the old checksum,
    /// `103a0436a5b1630ce94a6a8ee111f274`, exactly. Recorded in the Linux container (amd64).
    /// **Re-recorded 2026-10-22** when the default base-quality multiplier moved from 1.0 to 2.0:
    /// only the `##pvcParameters` line moved, carrying the parameters file's new checksum; with
    /// the old one put back it reproduces the old checksum, `7083dc4bede34f2e7f197c194ff8df41`, exactly. Recorded in the
    /// Linux container (amd64).
    /// **Re-recorded 2026-10-23** when `PARALOG_POST`'s header description became "under the paralog
    /// rate fitted around it" for the windowed duplication rate: only that `##INFO` line moved; putting
    /// the old description back reproduces the old checksum, `658380c537ec2cd5964a6774d2deafd7`, exactly. Recorded
    /// in the Linux container (amd64).
    const WITH_ERRORS_CALLS_MD5: &str = "d717848f594200eb3b7c4e39ddeeb3ae";

    /// **The same cohort grown to 21 samples**
    /// ([`a_larger_varying_cohort_with_sequencing_errors_on_disk`]), the fit and the calls. **A
    /// different computation of the errors**: above 20 samples the fit keeps each sample's own block
    /// of the information rather than the whole matrix, and inverts it a block at a time — the path
    /// every cohort the caller is built for takes, and one a two-sample cohort never reaches. Its 22
    /// libraries' multipliers run 0.35 to 7.97 with errors of 1.66 to 3.01; its inbreeding
    /// coefficients carry none; its calls are the 22 records designed into it, one site a sample and
    /// the tract. Recorded in the Linux container (arm64, glibc) and the same on macOS (arm64).
    ///
    /// **Re-recorded 2026-10-11, twice,** for the two reasons [`FITTED_PARAMETERS_MD5`] gives.
    /// Recorded in the Linux container (arm64, glibc) only; it was `5e3dacbf558adc6787276897935dbbf8`, then `0cda09cc72a2a7bd9c46a4fe85388eff`.
    ///
    /// **Re-recorded 2026-10-16** for the four reasons [`FITTED_PARAMETERS_MD5`] gives. Here the
    /// tract's row has F 0.897 and a slip rate of 0.00032 (5.4 × 10⁻⁸ at C7, before the stutter
    /// prior), and the fallback F is that F. Recorded on macOS (arm64) only; it was
    /// `0573095bb8f1ed6608ed74d62051ade1`. **Again the same day** for the note text of plan step
    /// D8 alone, as [`FITTED_PARAMETERS_MD5`] says (the tract is again a class of 1 locus); it was
    /// `5abcd04153d5acc32463fcd751970024`. **A third time** for the conditional notes alone, as
    /// [`FITTED_PARAMETERS_MD5`] says; it was `33eab2535f49967c6459db376a745d01`.
    ///
    /// **Re-recorded 2026-10-16** (`snp_fit_convergence.md` plan step A2) for the note on
    /// `snp_indel_fit_starts` alone: its three starts ran to the limit before and still do, and
    /// putting the old note back gives the old checksum exactly. The calls do not move. Recorded on
    /// macOS (arm64) only; it was `e0902dab0b2d64082d42a40b134f6721`.
    ///
    /// **Re-recorded 2026-10-16 again** (plan step B1), for the reason [`FITTED_PARAMETERS_MD5`]
    /// gives: every start still at the limit (199, 199, 198 passes), with 4, 18 and 22 Newton steps
    /// kept and 17, 17 and 20 judging passes whose every try was refused. 26 lines move: two starts'
    /// rows, all 22 libraries' multipliers (now 0.057 to 8.85 where they were 0.35 to 7.97), and the
    /// ordinary-site prior's `reference_concentration`, 44.10 to 50.00, and
    /// `alternative_concentration_total`, 0.0416 to 0.0580. Recorded on macOS (arm64) only; it was
    /// `66397f1ecd8c0effcf46bb53e97126bc`.
    ///
    /// **Re-recorded 2026-10-17** (plan step B2), for the reason [`FITTED_PARAMETERS_MD5`] gives:
    /// every start still at the limit (199, 198, 198 passes), for the two reasons
    /// [`WITH_ERRORS_FITTED_PARAMETERS_MD5`]'s note gives. The best log-likelihood rises from
    /// −5,108.54 to −5,053.47. All 22 multipliers move, now 0.0032 (three libraries on their rate's
    /// lower bound) to 13.99 with errors of 1.31 to 14.76, where they were 0.057 to 8.85; the
    /// ordinary-site prior's `alternative_concentration_total` 0.0580 to 0.0590 (the reference one
    /// stays at about 50.00); one start's row; and in the tenth significant figure the
    /// repeat-tract row and fallback F. Recorded on macOS (arm64) only; it was
    /// `8451e9ebed2660b2f35549d4844875d0`.
    ///
    /// **Re-recorded 2026-10-17 again** (plan step B3), for the reason
    /// [`WITH_ERRORS_FITTED_PARAMETERS_MD5`] gives. The swap had fired thousands of times here (6,126
    /// in the review's count). Every start still runs to the limit (198, 199, 199), but the best
    /// log-likelihood rises from −5,053.47 to −5,045.64, the review's 200-pass figure for the
    /// corrected fit; 12 of the 22 libraries end with the mismapped rate below the ordinary one, and
    /// the "mismapped" class holds about 45% of the positions (the review's measurement of this
    /// fix) — on this fixture the two classes
    /// split ordinary positions into two error levels. All 22 multipliers move, now 0.0032 to 15.59
    /// with errors of 1.03 to 4.73; the 21 inbreeding coefficients 0 to 7.8 × 10⁻¹¹; the
    /// ordinary-site prior's concentrations 50.00 to 49.91 and 0.0590 to 0.0498; the starts' rows;
    /// and the repeat-tract row and fallback F in their eleventh significant figure. Recorded on
    /// macOS (arm64) only; it was `291463840e49dca2f779491725538af3`.
    ///
    /// **Re-recorded 2026-10-17 a third time** (plan step B4), for the reason
    /// [`FITTED_PARAMETERS_MD5`] gives: **every start converges**, in 84, 38 and 44 passes, where all
    /// three ran to the limit; the best log-likelihood is −5,045.753 against −5,045.638 after 198
    /// passes. Its errors come from the blocks, so the decrement over the determined parameters is
    /// solved through the arrow here. All 22 multipliers move, now 0.0032 to 15.52 with errors of
    /// 0.83 to 9.14; the 21 inbreeding coefficients 7.8 × 10⁻¹¹ to 3.2 × 10⁻¹⁰; the ordinary-site
    /// prior's concentrations 49.91 to 49.52 and 0.0498 to 0.0526; the starts' rows; and the
    /// repeat-tract row and fallback F in their eleventh significant figure. Recorded on macOS
    /// (arm64) only; it was `ecd5c80409b9b372d50015d6c765c297`.
    ///
    /// **Re-recorded 2026-10-17 again** (plan step D2), for the reason [`FITTED_PARAMETERS_MD5`]
    /// gives: the version line and the two informational tables, both pairs marked — the density's
    /// first shape at 0.063 against an invariant share of 0.153, the duplicated share on its bound
    /// of 0.05 — and only the fixed non-reference share carrying an error. The tables deleted and
    /// the version set back give the old checksum; the calls' did not move. Recorded on macOS
    /// (arm64) only; it was `c2d56c4117a633b23c475ee7e422a53b`.
    ///
    /// **Re-recorded 2026-10-17 once more**, when the Newton step's box solve began treating a
    /// parameter it has not moved yet, sitting within the golden section's stopping width of the
    /// end its step heads to, as on that end (review of 2026-10-17, finding 4). The solve reaches
    /// the same maximum by another path, so its steps differ in their last bits, and the points the
    /// fit tries them at with them: the starts' rows do not move (the same endings and passes), the
    /// fitted numbers move by at most 6 × 10⁻⁶ of themselves and 4 × 10⁻⁶ of their errors, the
    /// errors by at most 4.4 × 10⁻⁵ of themselves. With that rule switched off the old checksum
    /// comes back; the calls' did not move. Recorded on macOS (arm64) only; it was
    /// `728256994c0c90109ae7f4a15eeb00a9`.
    ///
    /// **Re-recorded 2026-10-17 for the Newton step's target**, when the Newton step stopped putting a number the
    /// golden-section search maximises on an end of its interval: its target now stops half the
    /// search's stopping width short of each end (`Coordinate::newton_range` in `joint::fit`), the
    /// fix for two nearly fully inbred tomato samples whose excess a step put on exactly 1. Here
    /// the libraries' mismapped-position rates sit at their floor of 10⁻⁴, 3.8 × 10⁻¹⁰ above it,
    /// and the steps heading to it were clamped onto 10⁻⁴ where they now stop 5 × 10⁻¹⁰ above it;
    /// the tries along those steps differ from pass 40 of the first start. The starts' rows do not
    /// move (84, 38, 44 passes, all converged), nor the 21 inbreeding coefficients; the 22
    /// multipliers move by at most 0.00026 of their standard errors (4.5 × 10⁻⁴ of themselves),
    /// their errors by at most 0.2%; the ordinary-site prior's concentrations in their fifth
    /// figure; the informational tables' numbers by at most 0.09% (the duplicated class's share
    /// with carriers above the stated frequency). Recorded on macOS (arm64) only; it was
    /// `c548c669d12d20097217ff03e1299e3a`.
    ///
    /// **Re-recorded 2026-10-17** when the tract ladder's stated concentration went from one to ten
    /// (`doc/devel/reports/ng_str_single_sample_prior_2026-10-17.md`): this fit has a fitted median,
    /// so the number is not read, and only the repeat-tract note naming the built-in value moved
    /// (0.5 to 1/11). Recorded on macOS (arm64) only; it was `e0902dab0b2d64082d42a40b134f6721`.
    ///
    /// **Merged 2026-10-18 again**, as [`FITTED_PARAMETERS_MD5`]'s last note: main's two note
    /// passages alone; putting the old two back gives the branch's
    /// `abdbf13bab51eb6d15ba9b57944d214d` exactly. Recorded on the macOS host (arm64); it was
    /// `d623fef0663ce22a276fe551882dc06e` on main.
    ///
    /// **Re-recorded 2026-10-18** when the tract ladder's second rung got its own stated
    /// between-locus share (format 5, `doc/devel/reports/ng_str_single_sample_prior_2026-10-17.md`
    /// §8): this fixture's fit has a median of its own, so both keys carry it and no call moves.
    /// The file moves by `format_version`, the new key and the notes naming the two built-in
    /// values, and by nothing else. Recorded on the macOS host (arm64); it was `c47b81917edfbfa6f5c877eeae6087db`.
    /// **Re-recorded 2026-10-18** when the built-in repeat-tract outlier weight went from 0.20 to
    /// 0.07 (parameters file format 6): the file moves by `format_version` and the weight's line,
    /// and by nothing else. Recorded on the macOS host (arm64); it was `48480f4ee21a23f1580c059dec5e6bc5`.
    /// **Re-recorded 2026-10-22** when the default base-quality multiplier moved from 1.0 to 2.0:
    /// no number moved, because every read group here is fitted; only the file's prose about a
    /// `defaulted` multiplier did. Putting the old prose back reproduces the old checksum,
    /// `ced2576a45b958bd719a0e8720e7b9bb`, exactly. Recorded in the Linux container (amd64).
    /// **Re-recorded 2026-10-25** when `DECREMENT_THRESHOLD` went from 0.001 to 0.01
    /// (`fit_precision.md` §2, measured on 1, 8 and 62 tomato samples first): this fit's only start
    /// stops at 49 passes instead of 84, and its numbers and standard errors move within their
    /// errors — the largest, read group 15's multiplier, 0.47 → 0.58 against an error of 0.83.
    /// Setting the threshold back to 0.001 reproduces the old checksum,
    /// `b10618abc9bb00b0e089bfd64ae06858`, exactly. Recorded on the macOS host (arm64).
    const OF_21_SAMPLES_FITTED_PARAMETERS_MD5: &str = "aef60ff6d9816096e75a960d9ca5d0e8";
    /// **Re-recorded 2026-10-07** with the strand penalty of [`CALLS_MD5`]'s last note: the same 22
    /// records and genotypes, the SNPs' `QUAL` lower. It was `4a2c9b7ef87c94f09d55261bb290bf57`.
    /// Recorded in the Linux container (arm64, glibc) only.
    ///
    /// **Re-recorded 2026-10-10** when the hidden-duplication filter began scoring STR records from
    /// their read split as well as coverage (`hidden_paralog_filter.md` §3.2), measured on 63 tomato
    /// accessions first. The fit does not move. Only the STR record's `PARALOG_LR` moves (in the
    /// two-sample cohort, -1.11 to -19.24: its heterozygote's reads split 8 to 7, as a real
    /// heterozygote's do); no record, genotype or filter decision changes. It was `ec9aaa04a0d365fdaa468eb8eeda7c40`.
    /// Recorded in the Linux container (arm64, glibc) only.
    ///
    /// **Re-recorded 2026-10-10, again,** with the tract prior of [`CALLS_MD5`]'s last note: only the
    /// repeat tract's `QUAL` moves. It was `102696ee67202548f35a28a698dbf244`. Recorded in the
    /// Linux container (arm64, glibc) only.
    ///
    /// **Re-recorded 2026-10-11** when STR bundles began to be called (`str_bundle_loci.md` §7.5):
    /// the header's `RU` and `PERIOD` became `Number=.`, `STR`'s description names bundles, and
    /// `REPCNB` is declared. Those four header lines are all that moved — putting the old four back
    /// reproduces the old checksum exactly; the fixture holds no bundle, and no record moves. It
    /// was `fcc1a04573dbbdb08ac80127c554dc94`. Recorded in the Linux container (arm64, glibc) only.
    ///
    /// **Re-recorded 2026-10-12** with the partial-read change of [`WITH_ERRORS_CALLS_MD5`]'s
    /// last note: only the repeat tract's `QUAL` (37.2 to 37.3) and `AF` (0.023717 to 0.023718)
    /// move, and scoring every partial from the left again restores the old checksum.
    /// It was `fcc1a04573dbbdb08ac80127c554dc94`. Recorded in the Linux container (arm64,
    /// glibc) only.
    ///
    /// **Merged 2026-10-12**, as [`WITH_ERRORS_CALLS_MD5`]'s last note: the header is the bundle
    /// branch's and the repeat tract's record main's to the byte (`QUAL` 37.3, `AF=0.023718`). It
    /// was `08a675fd6ba4bc8a627fa2bd2118b415` on str-bundles and `efa3947cb4eacd394de0425ce87c40ed`
    /// on main. Recorded in the Linux container (arm64, glibc) only.
    ///
    /// **Re-recorded 2026-10-16** with the fit above (plan step D4): only the repeat tract at
    /// chrV:201 moves, QUAL 37.3 to 78.6, AF 0.023718 to 0.023811, the heterozygote's GQ 24 to 65
    /// and every other sample's GQ up by 2 or 3; no genotype. Recorded on macOS (arm64) only; it
    /// was `5c2414e47ecc6748f12db076d7513de1`.
    ///
    /// **Re-recorded 2026-10-16** with the fit above (`snp_fit_convergence.md` plan step B1; plan
    /// step A2 moved neither this nor the fit's numbers): 22 records to 23. One appears, a
    /// heterozygote at chrV:386 at QUAL 6.0 in one sample whose reads split 18 to 2; QUAL moves in
    /// 21 of the 22 that both hold, by −108.5 to +212.2, with the libraries' multipliers. No
    /// genotype of those 22 moved. Recorded on macOS (arm64) only; it was
    /// `b23487bd69f9d39f1943fbe9350a7a3f`.
    ///
    /// **Re-recorded 2026-10-17** with the fit above (plan step B2): 23 records to 24, one added at
    /// chrV:396 (QUAL 31.2); QUAL moves in 21 of the 23 both hold, by −96.0 to +325.2, with the
    /// multipliers; one genotype of 483 changes, the last sample's at chrV:55, 0/0 to 0/1; GQ moves
    /// in 10 records. Recorded on macOS (arm64) only; it was `64ac9af7ecaf0b2c603984e9b4420e6f`.
    ///
    /// **Re-recorded 2026-10-17 again** with the fit above (plan step B3): 24 records still, two
    /// lost (chrV:386 at QUAL 0.0 and chrV:396 at 31.2) and two added (chrV:483 at 26.8 and chrV:590
    /// at 45.3); QUAL moves in 21 of the 22 both hold, by −365.0 to +363.8, with the multipliers;
    /// three genotypes of 462 change — sample 7 at chrV:285 and chrV:377, 0/0 to 0/1, and the last
    /// sample at chrV:55 back to 0/0, as before step B2 — and GQ moves in almost every record.
    /// Recorded on macOS (arm64) only; it was `bfaff38cf2b6c5783ba280f98a177bb8`.
    ///
    /// **Re-recorded 2026-10-17 a third time** with the fit above (plan step B4): 24 records to
    /// 23, chrV:483 (QUAL 26.8) lost; QUAL moves in all 23 others, by under 3.4 but for chrV:354,
    /// 516.8 to 318.4; most samples' GQ by one; no genotype of 483 changes. Recorded on macOS
    /// (arm64) only; it was `aba287e8691aeb67c507d15005379ee1`.
    ///
    /// **Re-recorded 2026-10-17 again** with the fit above (the Newton target kept off the ends of
    /// the golden-section numbers): only the AF of two records moves in its sixth decimal, chrV:285
    /// 0.039376 to 0.039375 and chrV:377 0.048299 to 0.048298; no QUAL, genotype or GQ. Recorded on
    /// macOS (arm64) only; it was `c3e72793116b01d57ca63cbc3bbd2c16`.
    ///
    /// **Re-recorded 2026-10-16** when the hidden-duplication filter's VCF id became `hiddenDup`
    /// (`doc/devel/ng/spec/command_line.md` §4.3): its `##FILTER` header line is all that moved —
    /// putting `hiddenParalog` back reproduces the old checksum exactly. It was `b23487bd69f9d39f1943fbe9350a7a3f`, the value `main` recorded for the fit above.
    /// Recorded on the macOS host (arm64).
    ///
    /// **Re-recorded 2026-10-16** when a calling run stopped writing a parameters file beside its
    /// VCF (`doc/devel/ng/spec/command_line.md` §7): `##source` reads `pvc` where it read `ng`, and
    /// `##parametersFile=calls.parameters.toml` became `##pvcLog` and `##pvcParameters`. Those
    /// header lines are all that moved — putting the old two back reproduces the old checksum
    /// exactly. It was `45a981f5c11fde06848374dcdfe574a3`. Recorded on the macOS host (arm64).
    ///
    /// **Re-recorded 2026-10-17** when the coverage histogram's depth axis came to be fitted from
    /// every window through a fixed relative-depth grid (`doc/devel/ng/spec/window_coverage.md`
    /// §3.4): `PARALOG_LR` and `PARALOG_POST` are all that moved — with both taken out, the old
    /// and new files are identical line for line, the filter's header line included, and all
    /// 22 records are kept by both. The fixture is under 10,000 windows, so the old axis was
    /// already fitted from every window, and the move is the grid's rounding of the median and the
    /// re-cut alone. It was `23fb1a94c44118ecf249a68afaa67979`. Recorded on the macOS host (arm64).
    ///
    /// **Merged 2026-10-18**, main into `snp-fit-convergence`, as [`CALLS_MD5`]'s last note: only
    /// main's header lines, filter id and `PARALOG_LR`/`PARALOG_POST` differ from the branch's VCF;
    /// the branch's 23 records stay, with their genotypes, QUAL and GQ. It was
    /// `8c2025eaadcaf8f5e00a04688e2bd668` on the branch and `9d3c1c6c8243c41bb251d0d6cf7de7bb` on
    /// main. Recorded on the macOS host (arm64).
    ///
    /// **Re-recorded 2026-10-18**, merging the tract ladder's stated concentration (one to ten,
    /// `doc/devel/reports/ng_str_single_sample_prior_2026-10-17.md`) onto the above: only the
    /// `##pvcParameters` md5 moved — putting the old one back reproduces `9d3c1c6c8243c41bb251d0d6cf7de7bb` exactly.
    /// Recorded on the macOS host (arm64).
    ///
    /// **Merged 2026-10-18 again**, as [`CALLS_MD5`]'s last note: the `##pvcParameters` md5 alone;
    /// putting the branch's back gives the branch's `11ed05df550c9f1038bb07edf4975adc` exactly. It
    /// was `61cc1dbcd0bd87a01eadbdfadd64bd4e` on main. Recorded on the macOS host (arm64).
    ///
    /// **Re-recorded 2026-10-18** with the fit above: only the `##pvcParameters` md5 moved —
    /// putting the old one back reproduces the old checksum exactly. It was `2b753fd7c6d6e453645bb4b31fe6a364`. Recorded on
    /// the macOS host (arm64).
    /// **Re-recorded 2026-10-18** with the outlier weight above: besides the `##pvcParameters` md5,
    /// the one STR record (`chrV:201`) moves in QUAL and in its heterozygous sample's GQ — more
    /// confident at the lower floor — with every genotype unchanged. Recorded on the macOS host
    /// (arm64); it was `2883bcaad4b3f7552c384bf5fcf36335`.
    /// **Re-recorded 2026-10-21** when the hidden-duplication filter's allele term became a
    /// beta-binomial (`DEFAULT_ALLELE_DISPERSION`, `doc/devel/reports/ng_dup_filter_allele_balance_2026-10-21.md`):
    /// only all fifteen records' `PARALOG_LR`, each by 0.03 to 2.4 toward zero (chrV:201 −18.78 → −16.73 the most) moved; every record still passes, and `PARALOG_POST`, the header and every
    /// other column are unchanged. Setting the dispersion to zero reproduces `65c158aec014804780d5560220169f13` exactly.
    /// Recorded on the macOS host (arm64).
    /// **Re-recorded 2026-10-22** when the repeat-tract periodicity check was removed
    /// (`doc/devel/ng/spec/candidate_alleles_ssr.md` §7): only the header's
    /// `##FILTER=<ID=notPeriodic,…>` line went; putting it back reproduces the old checksum,
    /// `39efa7a9f70a5a07eea563e2f67df0c4`, exactly. Recorded in the Linux container (amd64).
    /// **Re-recorded 2026-10-22** when the default base-quality multiplier moved from 1.0 to 2.0:
    /// only the `##pvcParameters` line moved, carrying the parameters file's new checksum; with
    /// the old one put back it reproduces the old checksum, `2c68478585fa63e9b890cc30b3f33f9c`, exactly. Recorded in the
    /// Linux container (amd64).
    /// **Re-recorded 2026-10-23** when `PARALOG_POST`'s header description became "under the paralog
    /// rate fitted around it" for the windowed duplication rate: only that `##INFO` line moved; putting
    /// the old description back reproduces the old checksum, `bcd2e92c326cef11be6bcfb98798891b`, exactly. Recorded
    /// in the Linux container (amd64).
    /// **Re-recorded 2026-10-25** with [`OF_21_SAMPLES_FITTED_PARAMETERS_MD5`]'s last note: every
    /// record and every genotype is unchanged; the records' `QUAL` (chrV:32 220.1 → 195.7 the most,
    /// the rest by 1 or 2), `AF` in its fifth decimal, and the `##pvcParameters` line move. Setting
    /// the threshold back reproduces `513fb4a876d6f1d82b36a77171a4d825` exactly. Recorded on the
    /// macOS host (arm64).
    const OF_21_SAMPLES_CALLS_MD5: &str = "ae8b4f7101b1a69eb2ec1972e57aa69c";

    /// How many samples [`OF_21_SAMPLES_FITTED_PARAMETERS_MD5`] is of: one past the largest cohort
    /// whose errors come from the whole matrix (`fit::information::FULL_MATRIX_SAMPLES`, 20).
    const SAMPLES_PAST_THE_WHOLE_MATRIX: usize = 21;

    /// The 21-sample cohort.
    fn a_cohort_past_the_whole_matrix() -> AVaryingCohort {
        a_larger_varying_cohort_with_sequencing_errors_on_disk(SAMPLES_PAST_THE_WHOLE_MATRIX - 2)
    }

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
        let (psps, parameters) = walked_and_fitted(cohort);
        let calls = cohort.directory.path().join("calls.vcf");
        run_call_from_psps(&CallFromPspsArgs {
            reference: cohort.reference.clone(),
            catalog: Some(cohort.catalog.clone()),
            psps: cohort
                .samples
                .iter()
                .map(|sample| psp_path_for(&psps, sample))
                .collect(),
            output: calls.clone(),
            force: false,
            log: None,
            tmp_dir: None,
            parameters: Some(parameters.clone()),
            defaults: false,
            ploidy: None,
            max_locus_span: DEFAULT_MAX_COHORT_LOCUS_SPAN,
            max_generic_alleles: DEFAULT_MAX_CANDIDATE_ALLELES.get(),
            // The shipped false-discovery target, so the hidden-duplication filter scores and
            // writes its two fields: it computes with `ln`, `exp` and `ln_1p` as well.
            hidden_dup_fdr: 0.01,
            // Off: this test pins that platforms write the same bytes, over a fixture whose three
            // variant sites include low-quality ones; the threshold is pinned elsewhere.
            min_site_quality: 0.0,
            hidden_dup_tag: false,
            hidden_dup_window: 0,
            locus_builder_regions_len: None,
            psp_prefetch_bytes: crate::run::psp_prefetch::DEFAULT_PSP_PREFETCH_BUDGET_BYTES,
            threads: 0,
            max_reads_per_position: None,
            max_depth: crate::run::depth_ceiling::DEFAULT_MAX_READ_GROUP_DEPTH,
            explain_loci: None,
        })
        .expect("the cohort calls with its own fit");
        (parameters, calls)
    }

    /// Walk the cohort and fit it. Hands back the psps' directory and the parameters file's path.
    fn walked_and_fitted(cohort: &AVaryingCohort) -> (PathBuf, PathBuf) {
        let directory = cohort.directory.path();
        let psps = directory.join("psps");
        run_generate_psps(&GeneratePspsArgs {
            call_satellites: false,
            generic_path_only: false,
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
            max_str_bundle_len: DEFAULT_MAX_STR_BUNDLE_LEN,
            min_purity: DEFAULT_MIN_PURITY,
            max_depth: crate::locus_generation::pileup::DEFAULT_MAX_READS_PER_POSITION,
        })
        .expect("the cohort walks into psps");

        let parameters = directory.join("fitted.parameters.toml");
        run_estimate_parameters(&EstimateParametersArgs {
            reference: cohort.reference.clone(),
            catalog: Some(cohort.catalog.clone()),
            psps: vec![psps.clone()],
            output: parameters.clone(),
            force: false,
            log: None,
            ploidy: 2,
            inbreeding: None,
            str_param_estimates_at_once: None,
            skip_contamination: false,
            threads: 0,
        })
        .expect("the cohort fits");
        (psps, parameters)
    }

    /// One whole run — a fresh cohort walked, fitted and called — inside a pool of `threads`
    /// threads. Hands back the parameters file's text and the comparable VCF.
    fn one_run_at(threads: usize) -> (String, String) {
        one_run_of(a_varying_cohort_on_disk, threads)
    }

    /// The same, over the cohort `cohort_of` builds.
    fn one_run_of(cohort_of: fn() -> AVaryingCohort, threads: usize) -> (String, String) {
        let cohort = cohort_of();
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

    /// **The fit's standard errors have their recorded checksums** — recorded on arm64 Linux (glibc)
    /// and arm64 macOS; x86_64 not yet compared: the same cohort with sequencing errors in its
    /// reads, whose file carries a standard error on every library's multiplier.
    #[test]
    fn a_fit_with_standard_errors_writes_the_same_bytes_on_every_platform() {
        let (fitted, called) = one_run_of(a_varying_cohort_with_sequencing_errors_on_disk, THREADS);
        let multipliers_with_an_error = fitted
            .lines()
            .filter(|line| line.contains("error_probability_multiplier = { value = "))
            .filter(|line| line.contains("standard_error = "))
            .count();
        assert_eq!(
            multipliers_with_an_error, 3,
            "each of the three libraries' multipliers carries the error this test exists to pin:\n\
             {fitted}"
        );
        let fitted_md5 = md5_hex(fitted.as_bytes());
        let calls_md5 = md5_hex(called.as_bytes());
        assert_eq!(
            (fitted_md5.as_str(), calls_md5.as_str()),
            (WITH_ERRORS_FITTED_PARAMETERS_MD5, WITH_ERRORS_CALLS_MD5),
            "the checksums moved; read this module's documentation before re-recording them.\n\
             --- fitted parameters ---\n{fitted}\n--- calls ---\n{called}",
        );
    }

    /// **Every start of the SNP/indel fit on the two-sample cohort with sequencing errors converges,
    /// in fewer than 50 passes** (`fit_precision.md` §2, amendment of 2026-10-17). On this cohort
    /// the data leave most parameters without a standard error — 12 of 16, as the review of
    /// 2026-10-17 counted — so the decrement a start stops on counts only the 4 that have one, the
    /// others profiled out. The three starts converge after 22, 13 and 13 passes; counting every
    /// parameter, all three ran to the 200-pass limit. The checksums above were this rule's only
    /// guard before.
    #[test]
    fn every_start_on_the_cohort_with_sequencing_errors_converges() {
        use crate::calling::parameters_file::{ParametersFile, StartOutcome};
        let cohort = a_varying_cohort_with_sequencing_errors_on_disk();
        let (_, parameters) = walked_and_fitted(&cohort);
        let text = std::fs::read_to_string(&parameters).expect("the fit wrote its file");
        let file = ParametersFile::from_toml(&text).expect("the fit's file reads back");
        let starts = file
            .fitted_from
            .snp_indel_fit_starts
            .expect("the fit lists its starts");
        assert_eq!(starts.len(), 3, "{starts:?}");
        for start in &starts {
            assert_eq!(start.ended, StartOutcome::Converged, "{starts:?}");
            assert!(start.passes < 50, "{starts:?}");
        }
    }

    /// **The errors computed a sample's block at a time have their recorded checksums too**: 21
    /// samples, past the 20 at which the fit stops keeping the whole matrix. Every library's
    /// multiplier carries an error.
    #[test]
    fn errors_computed_a_block_at_a_time_write_the_same_bytes_on_every_platform() {
        let (fitted, called) = one_run_of(a_cohort_past_the_whole_matrix, THREADS);
        let read_groups = SAMPLES_PAST_THE_WHOLE_MATRIX + 1;
        let multipliers_with_an_error = fitted
            .lines()
            .filter(|line| line.contains("error_probability_multiplier = { value = "))
            .filter(|line| line.contains("standard_error = "))
            .count();
        assert_eq!(
            multipliers_with_an_error, read_groups,
            "each of the {read_groups} libraries' multipliers carries an error:\n{fitted}"
        );
        let fitted_md5 = md5_hex(fitted.as_bytes());
        let calls_md5 = md5_hex(called.as_bytes());
        assert_eq!(
            (fitted_md5.as_str(), calls_md5.as_str()),
            (OF_21_SAMPLES_FITTED_PARAMETERS_MD5, OF_21_SAMPLES_CALLS_MD5),
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
        // **Every cohort**: each guards against a chunk size taken from the pool's width, and the
        // two with sequencing errors carry it through the standard errors, by the whole matrix and
        // a block at a time. A join order chosen by the pool cannot show at the three chunks these
        // censuses make; `fit`'s `a_census_of_many_chunks_fits_to_the_same_bits_at_any_pool_width`
        // holds that.
        for cohort_of in [
            a_varying_cohort_on_disk as fn() -> AVaryingCohort,
            a_varying_cohort_with_sequencing_errors_on_disk,
            a_cohort_past_the_whole_matrix,
        ] {
            let at_four = one_run_of(cohort_of, THREADS);
            for threads in [1, 7] {
                let other = one_run_of(cohort_of, threads);
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
}
