use super::*;

/// **The same tract with the slippage round switched on** — `rounds` re-fit rounds at the
/// given pull-backs, everything else the shipped configuration. What the re-fit fixtures
/// call, and nothing else uses.
fn call_tract_with_refit(
    observations_of_each_sample: &[&[SequenceObservation]],
    rounds: u32,
    direction_and_fall_off_pull_back: f64,
    level_pull_back: f64,
    scratch: &mut CallingScratch<StutterSubstitutionScratch>,
) -> LocusInference {
    let detail = tract_detail();
    let per_sample: Vec<SsrSampleEvidence<'_>> = observations_of_each_sample
        .iter()
        .map(|observations| SsrSampleEvidence::new(observations, &detail, TEST_TRACT_LEN))
        .collect();
    let repeat_counts = tract_repeat_counts();
    let evidence = LocusEvidence::ssr(locus_region(), &per_sample, &detail, &repeat_counts);
    let inbreeding = vec![outbred(); per_sample.len()];
    let strata = tract_strata();
    let substitution = tract_substitution_rates();
    let calibration = tract_libraries();
    let parameters =
        uncontaminated_run(&calibration, &inbreeding, &strata, &substitution, diploid());
    let config = CallingLoopConfig {
        slippage_refit: crate::calling::inference::SlippageRefitConfig {
            max_rounds: rounds,
            direction_and_fall_off_pull_back_pseudocounts: direction_and_fall_off_pull_back,
            level_pull_back_slipped_reads: level_pull_back,
            round_convergence_threshold:
                crate::calling::inference::DEFAULT_ROUND_CONVERGENCE_THRESHOLD,
        },
        ..CallingLoopConfig::DEFAULT
    }
    .validate()
    .expect("a built setting with legal pull-backs");
    shipped_arm().call_locus(&evidence, &parameters, tract_alleles(), &config, scratch)
}

/// The three numbers the tract was last scored under at `(read group 0, candidate 0)` —
/// the one cell every re-fit fixture's reads land in.
fn effective_cell_0_0(
    scratch: &CallingScratch<StutterSubstitutionScratch>,
) -> crate::parameter_estimation::joint::ssr_fit::Slippage {
    scratch
        .tract_fits()
        .effective_slippage_of_cell(ReadGroupId(0), 0)
        .expect("the fixture's fit reaches this cell")
}

/// **A tract with planted heavy slippage pulls its level toward the reads' own rate and
/// is held short of it by the 20-slipped-reads pull-back**, with the expected numbers
/// computed by hand (spec §5.1).
///
/// One sample: 30 reads at the reference's six repeats and 10 reads at five — one whole
/// unit short, a length no candidate carries, so under the near-certain `0/0` call every
/// one of them is a one-unit contraction slip. The frozen cell (read group 0, stratum 6)
/// says `level 0.04, shorter_share 0.83, fall_off 0.25`.
///
/// The hand arithmetic, with the genotype posterior ≈ 1 on `0/0` (thirty clean reads
/// against ten one-step slips leave the alternatives tens of nats behind, so every weight
/// below is within ~10⁻⁶ of the whole read count):
///
/// - attributed weight ≈ 40, slipped ≈ 10, expected slips ≈ 40 × 0.04 = 1.6;
/// - level multiplier = (10 + 20) / (1.6 + 20) = 1.3888…, so the level moves
///   0.04 → **0.0555…** — toward the reads' own rate of 10/40 = 0.25, and nowhere near
///   it: the pull-back holds it at just over a fifth of the way;
/// - direction split = (10 + 50 × 0.83) / (10 + 50) = **0.8583…**;
/// - every slip is one unit, so the raw fall-off estimate is 0 and the blend is
///   (0 + 50 × 0.25) / 60 = **0.2083…**.
///
/// The second round attributes the same reads under the re-fitted table, lands on the
/// same pooled counts — the expected-slips denominator reads the **frozen** level, so the
/// multiplier does not compound — and stops for moving less than the threshold. Which is
/// also the early stop of the rounds' own rule: one round adopted of the three allowed.
#[test]
fn planted_slippage_moves_the_level_toward_the_reads_but_the_pull_back_holds_it() {
    let sample = [
        tract_reads(TRACT_CANDIDATE_REPEATS[0], 30),
        tract_reads(TRACT_CANDIDATE_REPEATS[0] - 1, 10),
    ];
    let mut scratch = worker_scratch();
    let inference = call_tract_with_refit(&[&sample], 3, 50.0, 20.0, &mut scratch);

    assert_eq!(called_alleles(&inference, 0), vec![0, 0]);
    let refitted = effective_cell_0_0(&scratch);
    assert!(
        (refitted.level - 0.04 * (30.0 / 21.6)).abs() < 1e-3,
        "level {} against the hand-computed 0.05555…",
        refitted.level
    );
    assert!(
        refitted.level < 0.25,
        "the pull-back must hold the level short of the reads' own rate"
    );
    assert!((refitted.shorter_share - 51.5 / 60.0).abs() < 1e-3);
    assert!((refitted.fall_off - 12.5 / 60.0).abs() < 1e-3);

    // The rounds stopped early — one adopted of the three allowed — and the one adopted
    // round is visible as the one extra table build.
    assert_eq!(scratch.slippage_refit_rounds(), 1);
    assert_eq!(scratch.emission_cost().emission_builds, 2);
}

/// **Zero pull-back is the free setting: the numbers go to the reads' own rate**
/// (HipSTR's behaviour, spec §5.1). Same reads as the planted fixture, so the hand
/// arithmetic differs only in the pull-backs: the multiplier is 10 / 1.6 = 6.25, the
/// level 0.04 × 6.25 = **0.25** — exactly 10 slipped reads in 40 — the direction split
/// **1.0** (every slip was a contraction) and the fall-off **0** (every slip one unit).
#[test]
fn zero_pull_back_hands_the_numbers_to_the_reads() {
    let sample = [
        tract_reads(TRACT_CANDIDATE_REPEATS[0], 30),
        tract_reads(TRACT_CANDIDATE_REPEATS[0] - 1, 10),
    ];
    let mut scratch = worker_scratch();
    let inference = call_tract_with_refit(&[&sample], 3, 0.0, 0.0, &mut scratch);

    assert_eq!(called_alleles(&inference, 0), vec![0, 0]);
    let refitted = effective_cell_0_0(&scratch);
    assert!(
        (refitted.level - 0.25).abs() < 1e-3,
        "free level {} against the reads' own rate 0.25",
        refitted.level
    );
    assert!((refitted.shorter_share - 1.0).abs() < 1e-3);
    assert!(refitted.fall_off.abs() < 1e-3);
    assert_eq!(scratch.slippage_refit_rounds(), 1);
}

/// **A tract whose reads carry no slips stays at the stratum's values.** Four clean reads
/// expect 4 × 0.04 = 0.16 slips, so the level multiplier is 20 / 20.16 = 0.9921 — and the
/// shape has no slips to move at all. The largest move any cell's numbers make is the
/// multiplier over the largest frozen level the gather holds (the second slippage group's
/// seven-repeat cell, 0.105): 0.105 × 0.0079 ≈ 0.0008, **under the round threshold of
/// 10⁻³** — so the very first round converges, its numbers are discarded (production's
/// test-before-adopt), and the locus keeps the frozen values **exactly**: no round
/// adopted, no table rebuilt.
///
/// *Why four reads and not ten*: the stopping rule reads every cell, and at ten clean
/// reads the multiplier is 20 / 20.4 = 0.9804, which moves that same 0.105 cell by
/// 0.002 — over the threshold, so a round is adopted and the level eases off the frozen
/// value. That is the behaviour production shares (its multiplier moves 0.0196 against
/// the same 10⁻³ tolerance), not a defect; this fixture sits below the threshold so that
/// *stays frozen* can be asserted bit for bit.
#[test]
fn a_tract_with_no_slips_stays_at_the_stratum_values() {
    let sample = [tract_reads(TRACT_CANDIDATE_REPEATS[0], 4)];
    let mut scratch = worker_scratch();
    let inference = call_tract_with_refit(&[&sample], 3, 50.0, 20.0, &mut scratch);

    assert_eq!(called_alleles(&inference, 0), vec![0, 0]);
    let cell = effective_cell_0_0(&scratch);
    assert_eq!(cell.level, 0.04, "the frozen value, bit for bit");
    assert_eq!(cell.shorter_share, 0.83);
    assert_eq!(cell.fall_off, 0.25);
    assert_eq!(scratch.slippage_refit_rounds(), 0);
    assert_eq!(
        scratch.emission_cost().emission_builds,
        1,
        "a converging first round rebuilds nothing"
    );
}

/// **An ordinary site ignores the re-fit structurally** (spec §5.1's closing paragraph):
/// the same SNP called with three rounds asked for and with none is the same call, the
/// same single table build, and a round count of zero.
#[test]
fn an_ordinary_site_ignores_the_slippage_rounds() {
    let (alleles, _) = generic_locus(1);
    let carrier = [observation(1, 0, 8, 4, 4)];
    let per_sample = [called_sample(&carrier)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let libraries = [ReadGroupCalibration::unscaled()];
    let inbreeding = [outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let substitution = std::collections::BTreeMap::new();
    let parameters = FrozenParameters::uncontaminated(
        &libraries,
        &inbreeding,
        human_like_seed(),
        &strata,
        &substitution,
        diploid(),
    );
    let config = CallingLoopConfig {
        slippage_refit: crate::calling::inference::SlippageRefitConfig {
            max_rounds: 3,
            ..crate::calling::inference::SlippageRefitConfig::DEFAULT
        },
        ..CallingLoopConfig::DEFAULT
    }
    .validate()
    .expect("re-fit rounds are a built setting");
    let mut scratch = worker_scratch();
    let inference =
        shipped_arm().call_locus(&evidence, &parameters, alleles, &config, &mut scratch);

    assert_eq!(called_alleles(&inference, 0), vec![1, 1]);
    assert_eq!(scratch.slippage_refit_rounds(), 0);
    assert_eq!(scratch.emission_cost().emission_builds, 1);
}
