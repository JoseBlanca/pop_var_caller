//! **A cohort on disk that any of this binary's commands can be driven over** — a reference,
//! a catalog built on it, and two samples' indexed alignment files.
//!
//! **One copy, and all three commands use it.** `generate-psps` and `call-from-alignments` each
//! grew their own, and Milestone C of `doc/devel/ng/impl_plan/run_driver_psp_mode.md` recorded
//! the duplication as a debt; `call-from-psps` would have been the third, and it is the one
//! command that has to be driven over *the files another command wrote*, so it cannot have a
//! private cohort at all.
//!
//! **The two copies were not the same cohort**, which is why sharing them was worth doing
//! rather than merely tidy: `call-from-alignments`' gave both samples no reads, so the only
//! test that drove that whole command drove it over a cohort with nothing in it.
//!
//! **What the cohort is built to catch**: one sample carries reads and the other carries none,
//! so a walk that produced records and the analysed-but-empty case (spec §12.9) are both in
//! every run driven over it — and a wiring defect that walked the right sample over the wrong
//! ground cannot hide behind two empty files.

use std::path::PathBuf;

use tempfile::TempDir;

use crate::ng::repeat_catalog::StrRepeatCriteria;

/// A reference, a catalog built on it, and two samples' alignment files — with the temporary
/// directories that keep them alive.
///
/// **The directories are returned rather than kept**, because a `TempDir` deletes its tree when
/// it drops: a fixture that held them internally would hand back paths to files that no longer
/// exist.
pub(crate) struct ACohortOnDisk {
    /// The reference FASTA, and the directory holding it, the catalog and anything a run writes
    /// beside them.
    pub reference: PathBuf,
    /// The tandem-repeat catalog built on that reference.
    pub catalog: PathBuf,
    /// The samples' alignment files, in the order `zeta`, `alpha` — the first carrying three
    /// reads, the second none.
    pub alignments: Vec<PathBuf>,
    /// Where the reference and the catalog live; also the natural home for a run's output.
    pub directory: TempDir,
    /// Keeps `zeta`'s alignment file alive.
    pub zeta: TempDir,
    /// Keeps `alpha`'s alignment file alive.
    pub alpha: TempDir,
}

/// Build the cohort. Everything is written to fresh temporary directories.
///
/// # Panics
///
/// On any failure to build the fixture, which is a broken test rather than a finding.
pub(crate) fn a_cohort_on_disk() -> ACohortOnDisk {
    use crate::ng::read::input::test_fixtures::{
        FIXTURE_CONTIGS, header, indexed_named_bam, matching_contigs, read_group_for,
        read_named_with_length_in_read_group,
    };
    use crate::ng::reference_info::{ReferenceSource, read_reference_info_observing};
    use crate::ng::repeat_catalog::RepeatCatalogBuilder;
    use crate::ng::tandem_repeat::ScanParams;
    use crate::pileup::per_sample::cram_files::{ContigSpec, build_fasta};
    use noodles_sam::alignment::RecordBuf;

    let specs: Vec<ContigSpec> = FIXTURE_CONTIGS
        .iter()
        .map(|(name, length)| ContigSpec {
            name: (*name).to_string(),
            length: *length as u64,
        })
        .collect();
    let (directory, reference) = build_fasta(&specs).expect("a reference on disk");

    let catalog = directory.path().join("ref.fa.repeats.parquet");
    let criteria = StrRepeatCriteria::default();
    let mut builder = RepeatCatalogBuilder::create(
        &catalog,
        criteria,
        ScanParams {
            match_reward: 2,
            mismatch_penalty: 7,
            min_copies: 2,
        },
    )
    .expect("a catalog to build into");
    let reference_info = read_reference_info_observing(
        ReferenceSource::Fasta {
            fasta: reference.clone(),
            fai: None,
        },
        &mut builder,
    )
    .expect("the reference reads");
    builder
        .finish(&reference_info)
        .expect("the catalog is written");

    let with_sample = |sample: &str, file: &str, records: &[RecordBuf]| {
        indexed_named_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[(&read_group_for(sample), Some(sample))],
            ),
            records,
            file,
        )
    };
    let zeta_reads = [
        read_named_with_length_in_read_group("z-r0", 0, 5, 30, &read_group_for("zeta")),
        read_named_with_length_in_read_group("z-r1", 0, 20, 30, &read_group_for("zeta")),
        read_named_with_length_in_read_group("z-r2", 1, 40, 30, &read_group_for("zeta")),
    ];
    let (zeta, zeta_bam) = with_sample("zeta", "zeta.bam", &zeta_reads);
    let (alpha, alpha_bam) = with_sample("alpha", "alpha.bam", &[]);

    ACohortOnDisk {
        reference,
        catalog,
        alignments: vec![zeta_bam, alpha_bam],
        directory,
        zeta,
        alpha,
    }
}
