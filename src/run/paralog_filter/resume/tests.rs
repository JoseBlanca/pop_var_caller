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

/// **A manifest says whether the core floor was set on its own, and an older one is read as the
/// run it was** (code review 2026-11-24, Mi2): a resumed run given a new `--refine-window-floor`
/// moves the core floor with it only where the stopped run left it unset. A manifest written
/// before the core floor existed had none set; one written before this was recorded had it set
/// exactly where it differs from the margin floor.
#[test]
fn a_manifest_says_whether_the_core_floor_was_set_on_its_own() {
    let inbreeding = vec![InbreedingF::try_new(0.9).expect("a fraction"); 2];
    let coefficients: Vec<f64> = inbreeding.iter().map(|f| f.get()).collect();
    let ploidy = Ploidy::try_new(2).expect("two copies");
    let default = decide(NeighbourEvidenceSwitch::On, &coefficients, ploidy);
    let given = default.clone().with_core_floor(300_000);
    // Written, then read back with `drop` keys taken out of the manifest's JSON.
    let read_back = |decision: &NeighbourEvidenceDecision, drop: &[&str]| {
        let mut json = serde_json::to_value(RefinementParts::of(decision)).expect("to JSON");
        for key in drop {
            json.as_object_mut().expect("an object").remove(*key);
        }
        serde_json::from_value::<RefinementParts>(json)
            .expect("from JSON")
            .into_decision(&inbreeding)
            .expect("a decision")
    };
    for keys in [
        &[][..],
        &["core_floor_given"],
        &["core_floor_given", "core_floor_bases"],
    ] {
        let back = read_back(&default, keys);
        assert!(!back.core_floor_given, "{keys:?}");
        assert_eq!(
            back.with_window_floor(800_000).sizes.core_floor_bases,
            800_000
        );
    }
    for keys in [&[][..], &["core_floor_given"]] {
        let back = read_back(&given, keys);
        assert!(back.core_floor_given, "{keys:?}");
        assert_eq!(
            back.with_window_floor(800_000).sizes.core_floor_bases,
            300_000
        );
    }
    assert_eq!(read_back(&given, &[]), given);
}
