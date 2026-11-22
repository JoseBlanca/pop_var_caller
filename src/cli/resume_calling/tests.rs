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
        phase_min_path_agreement: None,
    })
    .expect_err("no such directory");
    assert!(matches!(
        failed,
        ResumeCallingCliError::Calling(CallingRunError::ScratchDirectoryNotFound { .. })
    ));
}
