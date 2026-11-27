use super::*;

// ────────────────────────────────────────────────────────────────────────────────
// D2 — the two cost invariants: what a locus pays for, and what a pass does not
// ────────────────────────────────────────────────────────────────────────────────

/// One locus's evidence over `per_sample_observations` samples, where entry *s* says how
/// many `(allele, read group)` rows sample *s* shows — **deliberately unequal**, because
/// samples do not have equal observation counts and a fixture built so they do is the one
/// shape that hides a three-way product (`spec/calling_em_loop.md` §13's test 5).
fn rows_of_each_sample(per_sample_observations: &[usize]) -> Vec<Vec<GenericObservation>> {
    per_sample_observations
        .iter()
        .enumerate()
        .map(|(sample, observations)| {
            (0..*observations)
                .map(|row| {
                    // Every row a different allele, and the read counts varied by sample so
                    // that no two rows are interchangeable. **What the fixture is for is the
                    // count, not the call** — and the calls it happens to produce run the
                    // other way from the obvious guess: the sample whose reads are spread
                    // over three alleles is called at 54.7 Phred and the one-row sample at
                    // 12.3, because spreading the reads is what makes a heterozygote's two
                    // alleles both visible.
                    observation(row as u16, 0, 4 + sample as u32, 2, 2)
                })
                .collect()
        })
        .collect()
}

/// **`candidates × Σ_s (observations in sample s)`, and `Σ_s` is the part that needs
/// saying** — asserted as the formula rather than as a literal, at two loci that take
/// different numbers of passes.
///
/// Three samples showing one, two and three observations over three candidate alleles:
/// the sum is `3 × (1 + 2 + 3)` = **18**. **The two wrong shapes this fixture separates,
/// both measured on it:** charging the first row's count for every row gives
/// `3 × 1 × 3` = **9**, and charging the locus's pooled total for every row gives
/// `3 × 6 × 3` = **54**. A fixture whose three samples showed three observations each
/// would report 27 under all three, which is exactly what the equal-count fixture beside
/// this one cannot tell apart.
///
/// **The same 18 at two passes and at four.** That is what "the table is paid for once"
/// means: a build hoisted inside the frequency loop returns identical genotypes and costs
/// the emission four times over, so nothing but the counter can tell the two apart.
#[test]
fn the_emission_count_is_candidates_times_the_sum_over_samples_at_every_pass_count() {
    let (alleles, _) = generic_locus(2);
    let rows = rows_of_each_sample(&[1, 2, 3]);
    let per_sample: Vec<GenericLocusSample<'_>> =
        rows.iter().map(|row| called_sample(row)).collect();
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled()];
    let inbreeding = vec![outbred(); 3];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let parameters = uncontaminated_run(
        &calibration,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    );

    let candidates = 3_u64;
    let observations: u64 = [1, 2, 3].iter().sum();
    let mut passes_seen = Vec::new();
    // **One scratch across both loci, which is the only shape a real worker has** — and
    // what makes the counter's own reset load-bearing. Measured with
    // `prepare_for_locus`'s reset deleted (2026-08-26): the second locus reports
    // `emission_builds: 2, emission_row_fills: 6, emission_evaluations: 36,
    // table_assemblies: 2, row_assemblies: 6` — every field doubled — and a fresh scratch
    // per locus hides it completely.
    let mut scratch = worker_scratch();
    for config in [capped_at(2), RunnableCallingLoopConfig::default()] {
        let inference = shipped_arm().call_locus(
            &evidence,
            &parameters,
            alleles.clone(),
            &config,
            &mut scratch,
        );
        passes_seen.push(inference.passes);
        assert_eq!(
            scratch.emission_cost(),
            EmissionCost {
                emission_builds: 1,
                emission_row_fills: 3,
                emission_evaluations: candidates * observations,
                table_assemblies: 1,
                row_assemblies: 3,
            },
            "one emission build and one fold over three rows, {candidates} candidates × \
             {observations} observations, at {} passes",
            inference.passes
        );
    }
    assert_eq!(
        passes_seen,
        vec![2, 4],
        "the two runs must differ in their pass count, or the invariant is untested"
    );
}

/// **No buffer of the worker's scratch moves or grows, however many passes the loop
/// takes** — the loop's zero-allocation invariant, as far as a crate that forbids
/// `unsafe` can observe it (`CallingScratch::buffer_fingerprints`).
///
/// The same locus is called twice on **one** scratch, capped at two passes and then at
/// the default, and every buffer comes back at the same address with the same length. A
/// `Vec` that grew during a pass would have moved its bytes; a per-pass buffer added to
/// the scratch would have changed a length.
///
/// **What this cannot see is a temporary allocated and dropped inside a pass**, which
/// leaves no trace anywhere — so that half is counted for real, in
/// `tests/ng_calling_loop_allocation.rs`. It lives in a test binary of its own because a
/// global allocator counts the whole process and the lib suite runs in parallel; this
/// test is the cheap guard that runs on every build.
#[test]
fn no_buffer_of_the_scratch_moves_or_grows_however_many_passes_the_loop_takes() {
    let (alleles, _) = generic_locus(2);
    let rows = rows_of_each_sample(&[1, 2, 3]);
    let per_sample: Vec<GenericLocusSample<'_>> =
        rows.iter().map(|row| called_sample(row)).collect();
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled()];
    let inbreeding = vec![outbred(); 3];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let parameters = uncontaminated_run(
        &calibration,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    );
    let mut scratch = worker_scratch();

    let two_passes = shipped_arm().call_locus(
        &evidence,
        &parameters,
        alleles.clone(),
        &capped_at(2),
        &mut scratch,
    );
    let after_two = scratch.buffer_fingerprints();

    let four_passes = shipped_arm().call_locus(
        &evidence,
        &parameters,
        alleles,
        &RunnableCallingLoopConfig::default(),
        &mut scratch,
    );
    let after_four = scratch.buffer_fingerprints();

    assert_eq!(
        (two_passes.passes, four_passes.passes),
        (2, 4),
        "the two runs must differ in their pass count, or the invariant is untested"
    );
    assert_eq!(
        after_two, after_four,
        "every buffer is where it was and as long as it was: a Vec that grew during a \
         pass would have moved its bytes"
    );
}

/// **And the same on the contaminated path, which is the one that does per-pass work
/// inside the loop.**
///
/// The test above calls a run with no fraction fitted, so its loop reads a table nothing
/// touches. With a fraction fitted the loop refills the batch copies, refills one sample's
/// contaminant frequencies per row and assembles the whole table again, at every pass —
/// three buffers written inside the loop and sized outside it, which is exactly what
/// `buffer_fingerprints` claims to measure and what nothing measured until this fixture.
#[test]
fn no_buffer_moves_or_grows_on_a_contaminated_locus_either() {
    let (alleles, _) = generic_locus(2);
    let rows = rows_of_each_sample(&[1, 2, 3]);
    let per_sample: Vec<GenericLocusSample<'_>> =
        rows.iter().map(|row| called_sample(row)).collect();
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled()];
    let fractions = [contaminated_at(0.05)];
    let inbreeding = vec![outbred(); 3];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let batching = SequencingBatches::all_together_over(1, 3);
    let parameters = contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);
    let mut scratch = worker_scratch();

    let fewer = shipped_arm().call_locus(
        &evidence,
        &parameters,
        alleles.clone(),
        &capped_at(2),
        &mut scratch,
    );
    let after_fewer = scratch.buffer_fingerprints();

    let more = shipped_arm().call_locus(
        &evidence,
        &parameters,
        alleles,
        &RunnableCallingLoopConfig::default(),
        &mut scratch,
    );
    let after_more = scratch.buffer_fingerprints();

    assert!(
        more.passes > fewer.passes,
        "the two runs must differ in their pass count, or the invariant is untested — \
         {} against {}",
        fewer.passes,
        more.passes
    );
    // **Sorted, because two of the buffers exchange pointers with every pass** — the M-step
    // swaps the cohort's expected copies with the previous pass's rather than copying
    // either, so the list's order depends on the parity of the pass count and its *contents*
    // are the claim. The uncontaminated fixture above compares unsorted only because both
    // its runs take an even number of passes.
    let mut after_fewer = after_fewer;
    let mut after_more = after_more;
    after_fewer.sort_unstable();
    after_more.sort_unstable();
    assert_eq!(
        after_fewer, after_more,
        "every buffer is where it was and as long as it was, over a loop that refills three \
         of them at every pass"
    );
}

/// **A wider locus does grow the buffers, and it is meant to** — the fingerprint is not
/// a claim that the scratch never allocates, only that a *pass* does not.
///
/// Without this, the test above passes against an implementation whose buffers never
/// change because nothing ever prepares them.
#[test]
fn a_wider_locus_than_the_worker_has_met_does_grow_the_scratch() {
    let (narrow, _) = generic_locus(1);
    let (wide, _) = generic_locus(3);
    let rows = rows_of_each_sample(&[1]);
    let per_sample: Vec<GenericLocusSample<'_>> =
        rows.iter().map(|row| called_sample(row)).collect();
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled()];
    let inbreeding = [outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let parameters = uncontaminated_run(
        &calibration,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    );
    let mut scratch = worker_scratch();

    let _ = shipped_arm().call_locus(
        &evidence,
        &parameters,
        narrow,
        &RunnableCallingLoopConfig::default(),
        &mut scratch,
    );
    let after_narrow = scratch.buffer_fingerprints();

    let _ = shipped_arm().call_locus(
        &evidence,
        &parameters,
        wide,
        &RunnableCallingLoopConfig::default(),
        &mut scratch,
    );
    let after_wide = scratch.buffer_fingerprints();

    assert_ne!(
        after_narrow, after_wide,
        "a locus of four alleles is wider than one of two, so the per-genotype buffers \
         have to grow — a fingerprint that never moved would mean nothing was prepared"
    );
}
