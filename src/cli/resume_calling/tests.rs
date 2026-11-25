use super::*;

/// **A directory where no run kept a resume point is refused**, naming the directory, before any
/// output is made.
#[test]
fn a_directory_with_no_resume_point_is_refused() {
    let directory = tempfile::tempdir().expect("a temporary directory");
    let output = directory.path().join("calls.vcf");
    let failed = run_resume_calling(&ResumeCallingArgs {
        scratch_dir: directory.path().to_path_buf(),
        output: output.clone(),
        threads: 0,
        force: false,
        log: Some(directory.path().join("calls.log")),
        keep_tmp_on_error: true,
        refine_window_bytes: None,
        refine_window_floor: None,
        refine_margin_density: None,
        refine_ignore_memory_limit: false,
        refine_core_floor: None,
        phase_min_path_agreement: None,
    })
    .expect_err("nothing to resume");
    assert!(matches!(failed, ResumeCallingCliError::NotResumable { .. }));
    assert!(!output.exists());
}

/// **A scratch directory that is not there is refused** as such.
#[test]
fn a_scratch_directory_that_is_not_there_is_refused() {
    let directory = tempfile::tempdir().expect("a temporary directory");
    let failed = run_resume_calling(&ResumeCallingArgs {
        scratch_dir: directory.path().join("gone"),
        output: directory.path().join("calls.vcf"),
        threads: 0,
        force: false,
        log: Some(directory.path().join("calls.log")),
        keep_tmp_on_error: false,
        refine_window_bytes: None,
        refine_window_floor: None,
        refine_margin_density: None,
        refine_ignore_memory_limit: false,
        refine_core_floor: None,
        phase_min_path_agreement: None,
    })
    .expect_err("no such directory");
    assert!(matches!(
        failed,
        ResumeCallingCliError::Calling(CallingRunError::ScratchDirectoryNotFound { .. })
    ));
}

/// The resume options, none of them given.
fn no_resume_options() -> ResumeCallingArgs {
    ResumeCallingArgs {
        scratch_dir: PathBuf::from("unused"),
        output: PathBuf::from("unused.vcf"),
        threads: 0,
        force: false,
        log: None,
        keep_tmp_on_error: false,
        refine_window_bytes: None,
        refine_window_floor: None,
        refine_margin_density: None,
        refine_ignore_memory_limit: false,
        refine_core_floor: None,
        phase_min_path_agreement: None,
    }
}

/// **A resumed run plans the windows a fresh run with the same options plans** (code review
/// 2026-11-24, Mi2). A fresh run's core floor is its `--refine-window-floor` unless
/// `--refine-core-floor` is given; a run resumed with a new window floor kept the stopped run's
/// core floor instead, so cores could shrink below the margins a fresh run would have stopped at.
#[test]
fn a_resumed_run_plans_the_windows_a_fresh_run_with_the_same_options_plans() {
    let inbreeding = vec![crate::types::InbreedingF::try_new(0.9).expect("a fraction"); 12];
    let fresh = |window_floor: u64, core_floor: Option<u64>| {
        calling_run::neighbour_evidence_decision(
            calling_run::RefineGenotypesArg::On,
            &inbreeding,
            crate::types::Ploidy::try_new(2).expect("two copies"),
            PathAgreement::DEFAULT,
            calling_run::RefineWindowArgs {
                window_bytes: crate::neighbour_evidence::windows::DEFAULT_WINDOW_BYTES,
                window_floor,
                core_floor,
                margin_density: calling_run::MarginDensity(None),
                window_max: None,
                core_max: None,
            },
        )
    };
    let resumed = |stopped, window_floor: Option<u64>, core_floor: Option<u64>| {
        the_decision_resumed(
            stopped,
            &ResumeCallingArgs {
                refine_window_floor: window_floor,
                refine_core_floor: core_floor,
                ..no_resume_options()
            },
        )
    };
    // The core floor never given: it follows the new window floor.
    assert_eq!(
        resumed(fresh(500_000, None), Some(800_000), None).sizes,
        fresh(800_000, None).sizes
    );
    // Given to the stopped run: it stays, as it would in a fresh run given both.
    assert_eq!(
        resumed(fresh(500_000, Some(300_000)), Some(800_000), None).sizes,
        fresh(800_000, Some(300_000)).sizes
    );
    // Given to this run.
    assert_eq!(
        resumed(fresh(500_000, None), Some(800_000), Some(400_000)).sizes,
        fresh(800_000, Some(400_000)).sizes
    );
    // Nothing given: the stopped run's windows.
    assert_eq!(
        resumed(fresh(600_000, None), None, None).sizes,
        fresh(600_000, None).sizes
    );
}
