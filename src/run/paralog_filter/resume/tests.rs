use super::*;
use crate::neighbour_evidence::stage::{NeighbourEvidenceSwitch, decide};
use crate::window_coverage::CoverageByGcHistogram;

/// A finished, empty spill in its own directory.
fn a_finished_spill(directory: &Path) -> SpillFile {
    let mut spill = SpillFile::in_directory(directory);
    spill.finish_writing().expect("an empty spill finishes");
    spill
}

fn the_metadata(inbreeding: &[InbreedingF]) -> VcfHeaderMetadata {
    let ploidy = Ploidy::try_new(2).expect("two copies");
    let coefficients: Vec<f64> = inbreeding.iter().map(|f| f.get()).collect();
    VcfHeaderMetadata::try_new(
        vec![
            HeaderContig {
                name: "chr1".to_string(),
                length: 1_000,
                md5: Some([7; 16]),
            },
            HeaderContig {
                name: "chr2".to_string(),
                length: 2_000,
                md5: None,
            },
        ],
        vec!["a".to_string(), "b".to_string()],
        "pvc call-from-psps --psp a.psp --psp b.psp".to_string(),
        "/genome.fa".to_string(),
        "pvc_calls.log".to_string(),
        "parameters.toml (md5 0123)".to_string(),
    )
    .expect("a valid header")
    .the_low_base_quality_test_ran()
    .the_neighbour_evidence_decision(
        decide(NeighbourEvidenceSwitch::On, &coefficients, ploidy)
            .with_window_bytes(123_456_789)
            .with_window_floor(250_000)
            .with_path_agreement(PathAgreement::try_new(0.7).expect("a share")),
    )
}

/// **A resume point comes back as the values calling handed the filter**, bit for bit: the
/// histograms, the inbreeding, what was asked of the filter, the header with its refinement
/// decision, the ploidy, and the spill's entry count.
#[test]
fn a_resume_point_is_read_back_as_it_was_written() {
    let directory = tempfile::tempdir().expect("a temporary directory");
    let spill = a_finished_spill(directory.path());
    let inbreeding = vec![
        InbreedingF::try_new(0.1 + 0.2).expect("a fraction"),
        InbreedingF::try_new(0.985_360_790_626_451_1).expect("a fraction"),
    ];
    let histograms = vec![
        SampleHistogram::Fitted(CoverageByGcHistogram {
            window_bp: 1_000,
            gc_bins: 2,
            depth_bin_width: 1.0 / 3.0,
            depth_bins: 1,
            windows_folded: 5,
            windows_under_the_floor: 1,
            depth_summed: 40,
            depth_rises: 2,
            counts: vec![1, 0, 3, 1],
        }),
        SampleHistogram::NoWindowFinalised,
    ];
    let asked_for = WhatTheOperatorAskedFor {
        target_fdr: TargetFdr::try_new(0.01).expect("a target"),
        tag_instead_of_dropping: true,
        window_bases: 50_000,
    };
    let metadata = the_metadata(&inbreeding);
    let ploidy = Ploidy::try_new(2).expect("two copies");

    assert!(!is_resumable(directory.path()));
    save(
        &spill,
        &histograms,
        &inbreeding,
        asked_for,
        &metadata,
        ploidy,
    )
    .expect("saved");
    assert!(is_resumable(directory.path()));

    let back = load(directory.path()).expect("read back");
    assert_eq!(back.histograms, histograms);
    assert_eq!(
        back.inbreeding
            .iter()
            .map(|f| f.get().to_bits())
            .collect::<Vec<_>>(),
        inbreeding
            .iter()
            .map(|f| f.get().to_bits())
            .collect::<Vec<_>>()
    );
    assert_eq!(back.asked_for, asked_for);
    assert_eq!(back.metadata, metadata);
    assert_eq!(back.ploidy, ploidy);
    assert_eq!(back.spill.entries_written(), spill.entries_written());
}

/// **A manifest from a build whose spill is laid out otherwise is refused**, rather than its
/// spill misread.
#[test]
fn a_resume_point_from_another_spill_layout_is_refused() {
    let directory = tempfile::tempdir().expect("a temporary directory");
    let spill = a_finished_spill(directory.path());
    let inbreeding = vec![
        InbreedingF::try_new(0.5).expect("a fraction"),
        InbreedingF::try_new(0.5).expect("a fraction"),
    ];
    save(
        &spill,
        &[
            SampleHistogram::NoWindowFinalised,
            SampleHistogram::NoWindowFinalised,
        ],
        &inbreeding,
        WhatTheOperatorAskedFor {
            target_fdr: TargetFdr::try_new(0.01).expect("a target"),
            tag_instead_of_dropping: false,
            window_bases: 0,
        },
        &the_metadata(&inbreeding),
        Ploidy::try_new(2).expect("two copies"),
    )
    .expect("saved");
    let path = directory.path().join(RESUME_FILE_NAME);
    let text = std::fs::read_to_string(&path).expect("the manifest");
    let changed = text.replace(
        &format!("\"spill_format\": {SPILL_FORMAT}"),
        &format!("\"spill_format\": {}", SPILL_FORMAT + 1),
    );
    assert_ne!(changed, text, "the manifest names the spill's layout");
    std::fs::write(&path, changed).expect("rewritten");

    assert!(matches!(
        load(directory.path()),
        Err(ResumeError::OtherFormat { what: "spill", .. })
    ));
}

/// **A directory where calling never ended holds no resume point**, and reading one says so.
#[test]
fn a_directory_without_a_manifest_is_not_resumable() {
    let directory = tempfile::tempdir().expect("a temporary directory");
    assert!(!is_resumable(directory.path()));
    assert!(matches!(
        load(directory.path()),
        Err(ResumeError::Read { .. })
    ));
}
