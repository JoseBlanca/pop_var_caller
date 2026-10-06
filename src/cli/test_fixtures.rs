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

use crate::repeat_catalog::StrRepeatCriteria;

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
    use crate::bam::cram_files::{ContigSpec, build_fasta};
    use crate::read::input::test_fixtures::{
        FIXTURE_CONTIGS, header, indexed_named_bam, matching_contigs, read_group_for,
        read_named_with_length_in_read_group,
    };
    use crate::reference_info::{ReferenceSource, read_reference_info_observing};
    use crate::repeat_catalog::RepeatCatalogBuilder;
    use crate::tandem_repeat::ScanParams;
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

/// The one contig the varying cohort is built on, and how long it is.
pub(crate) const VARYING_CONTIG: (&str, usize) = ("chrV", 600);

/// **Where the deliberate repeat tract sits**, 0-based and half-open, and what tiles it.
///
/// Ten copies of `GT` against a period-2 floor of six, with non-repetitive flanks either side,
/// so the segmentation types this stretch as a repeat tract and the rest of the contig as
/// ordinary sequence. **Without it the fixture exercises none of the record format's
/// repeat-tract half** — the motif, the two flanks, the per-allele repeat counts — and a psp
/// that stored `Generic` for every locus would pass the mode-equivalence oracle, which was
/// measured before this tract existed.
pub(crate) const TRACT: (usize, usize) = (200, 220);
/// The motif that tiles [`TRACT`].
pub(crate) const TRACT_MOTIF: &[u8; 2] = b"GT";

/// Where the cohort's samples carry a substitution, 0-based.
///
/// **One in each sample and neither in both**, so the two samples' stored files differ in what
/// they hold rather than only in the name in their headers — a defect that gave one sample's
/// observations to the other would otherwise be invisible, every record's two sample columns
/// being the same string.
pub(crate) const FIRST_SAMPLES_SUBSTITUTION: usize = 120;
/// The second sample's, which the first carries no read at.
pub(crate) const SECOND_SAMPLES_SUBSTITUTION: usize = 455;

/// A cohort whose samples **vary from the reference**, so that a run over it writes VCF records.
///
/// **[`a_cohort_on_disk`] cannot do this and that is not a fault in it.** Its reference is all
/// `A`s, so every base of it is one homopolymer run: the catalog routes the whole genome to the
/// repeat-tract path and a run over it writes no record at all. That is the right fixture for
/// asking what a command *writes beside* its VCF, and the wrong one for asking whether two
/// commands write the same VCF — two empty files are equal for the wrong reason.
///
/// **What this one is built to discriminate**, each closing a defect measured surviving the
/// oracle before it was added:
///
/// - **the two samples differ in what they carry**, not only in their names — one substitution
///   each, in different places, plus a tract length only the first sample varies at;
/// - **both strands are read and the alternative reads lean to one of them**, so a stored
///   forward-read count that was destroyed on write stops cancelling between the reference and
///   alternative reads;
/// - **the first sample declares two read groups and the second one**, so a walk numbers more
///   than one group from zero and the calling run's renumbering has something to do;
/// - **one locus is a repeat tract with a length variant in it**, so the record format's motif,
///   flanks and repeat counts are written and read rather than skipped.
pub(crate) struct AVaryingCohort {
    /// The reference FASTA.
    pub reference: PathBuf,
    /// The catalog built on it.
    pub catalog: PathBuf,
    /// The samples' alignment files, in the order of [`Self::samples`].
    pub alignments: Vec<PathBuf>,
    /// The samples' names, `one` and `two` and any further ones, in the order of
    /// [`Self::alignments`].
    pub samples: Vec<String>,
    /// Holds the reference, the catalog, and anywhere a run writes its output.
    pub directory: TempDir,
    /// **Keeps the samples' alignment files alive, and is read by nothing.** A `TempDir`
    /// deletes its tree when it drops, so a fixture that did not hand these back would return
    /// paths to files that no longer exist.
    pub _files: Vec<TempDir>,
}

/// **A reference sequence with no repeat in it**, so the segmentation routes it to the
/// SNP/indel path rather than to the repeat-tract path.
///
/// A fixed 32-bit congruential sequence, two bits a base, **rejecting any base equal to the one
/// before it** — which makes a homopolymer run of two impossible, let alone the eight copies
/// period 1 needs. Deterministic, so every run of the suite builds the same genome.
///
/// **What it does not guarantee is the longer periods**, whose floors are six copies at periods
/// 2 to 4, five at period 5 and four at period 6. Measured on the sequence this returns at
/// length 600: the longest tandem stretch is 5.5 copies of a period-2 motif, one copy short of
/// its floor, and nothing else reaches half a floor at any period. **That margin is one base
/// wide and nothing about the generator defends it**, which is why
/// `the_fixtures_ground_is_typed_as_its_doc_says` asserts the segmentation rather than trusting
/// this note: change the seed or the length and the assertion, not a reader, finds out.
fn a_reference_with_no_repeat_in_it(length: usize) -> Vec<u8> {
    const BASES: [u8; 4] = *b"ACGT";
    let mut bases = Vec::with_capacity(length);
    let mut state: u32 = 0x5EED_1234;
    while bases.len() < length {
        state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
        let base = BASES[((state >> 16) & 3) as usize];
        if bases.last() == Some(&base) {
            continue;
        }
        bases.push(base);
    }
    bases
}

/// **The reference the varying cohort is built on**: [`a_reference_with_no_repeat_in_it`] with
/// [`TRACT`] overwritten by ten copies of [`TRACT_MOTIF`].
///
/// Exposed so a test can read a base out of it and say what it expects a read carrying the
/// alternative to look like, rather than restating the sequence.
pub(crate) fn the_varying_cohorts_reference() -> Vec<u8> {
    let (name, length) = VARYING_CONTIG;
    let _ = name;
    let mut bases = a_reference_with_no_repeat_in_it(length);
    let (from, to) = TRACT;
    for (offset, base) in bases[from..to].iter_mut().enumerate() {
        *base = TRACT_MOTIF[offset % TRACT_MOTIF.len()];
    }
    bases
}

/// The base that is not `reference` — any of the other three will do, and picking
/// deterministically keeps the fixture reproducible.
///
/// **It says nothing about the neighbours**, so a read carrying this may hold a two-base
/// homopolymer the reference does not. That is harmless here because routing is computed from
/// the reference and never from the reads.
fn a_base_that_is_not(reference: u8) -> u8 {
    match reference {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    }
}

/// Build the varying cohort. Everything is written to fresh temporary directories.
///
/// **The stride is what makes the fixture call anything, and it was measured rather than
/// guessed.** At a stride of forty the same reads give 1.4 reads a locus and the run writes no
/// record at all — the cohort is analysed, the loci are built, and nothing about two reads is
/// evidence of a variant. At five it is twenty reads a position across the middle of the contig,
/// ten of them carrying the alternative, which is a heterozygote a defaults run scores.
///
/// # Panics
///
/// On any failure to build the fixture, which is a broken test rather than a finding.
pub(crate) fn a_varying_cohort_on_disk() -> AVaryingCohort {
    a_varying_cohort_built(false, 0)
}

/// **The same cohort with a sequencing error in about one read in four**, at an offset and to a base
/// hashed from the sample's name and the read, and never on the base that read's own sample varies
/// at (another sample's site, and the tract, can take one).
///
/// The cohort above has no error in any read, so the fit puts every library's error rate on its
/// lower bound and gives it no standard error, and its tract never mismatches. A file pinned from it
/// carries the format of the fit's errors and none of their arithmetic. **This one gives the fit
/// errors to compute**: each library's rate leaves its bound, so its multiplier carries a standard
/// error, and the tract's substitution count is no longer zero
/// (`cli::cross_platform_digests`, plan `fit_precision.md`, checkpoint E).
pub(crate) fn a_varying_cohort_with_sequencing_errors_on_disk() -> AVaryingCohort {
    a_varying_cohort_built(true, 0)
}

/// **The same, with `further_samples` samples beside the first two**, each read from one library
/// and carrying one substitution of its own, at a site away from the tract and from every other
/// sample's. For a cohort large enough that the fit computes its standard errors a sample's block
/// at a time rather than from the whole matrix, which it does above
/// `fit::information::FULL_MATRIX_SAMPLES` (20) samples.
pub(crate) fn a_larger_varying_cohort_with_sequencing_errors_on_disk(
    further_samples: usize,
) -> AVaryingCohort {
    a_varying_cohort_built(true, further_samples)
}

fn a_varying_cohort_built(with_sequencing_errors: bool, further_samples: usize) -> AVaryingCohort {
    use crate::read::input::test_fixtures::{
        FixtureReadGroup, header_with_read_groups, indexed_named_bam,
    };
    use crate::reference_info::{ReferenceSource, read_reference_info_observing};
    use crate::repeat_catalog::RepeatCatalogBuilder;
    use crate::tandem_repeat::ScanParams;
    use noodles_sam::alignment::record::MappingQuality;
    use noodles_sam::alignment::record::cigar::op::{Kind, Op};
    use noodles_sam::alignment::record::data::field::Tag;
    use noodles_sam::alignment::record_buf::data::field::Value;
    use noodles_sam::alignment::record_buf::{Cigar, QualityScores, Sequence};
    use noodles_sam::alignment::{RecordBuf, record::Flags};

    let (name, length) = VARYING_CONTIG;
    let bases = the_varying_cohorts_reference();

    let directory = tempfile::tempdir().expect("a temporary directory");
    let reference = directory.path().join("ref.fa");
    std::fs::write(
        &reference,
        format!(
            ">{name}\n{}\n",
            std::str::from_utf8(&bases).expect("ACGT is text")
        ),
    )
    .expect("the reference writes");

    let catalog = directory.path().join("ref.fa.repeats.parquet");
    let mut builder = RepeatCatalogBuilder::create(
        &catalog,
        StrRepeatCriteria::default(),
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

    // **Reads that match the reference except where the fixture says they vary**, so the read
    // filters keep them: a read of arbitrary bases would exceed the mismatch fraction and be
    // dropped, and the cohort would call nothing while every file looked full.
    let reads_of = |sample: &str,
                    substitution: usize,
                    shortens_the_tract: bool,
                    read_groups: &[String]|
     -> Vec<RecordBuf> {
        const READ: usize = 100;
        const STRIDE: usize = 5;
        let (tract_from, tract_to) = TRACT;
        (0..(length / STRIDE))
            .map(|index| {
                let start = index * STRIDE;
                let mut sequence = bases[start..(start + READ).min(length)].to_vec();
                // **Every second read carries the alternative**, which is a heterozygote — the
                // genotype a cohort of two can be called at.
                let carries_the_alternative = index % 2 == 0;
                if carries_the_alternative
                    && substitution >= start
                    && substitution < start + sequence.len()
                {
                    sequence[substitution - start] = a_base_that_is_not(bases[substitution]);
                }
                // **A tract shorter by one copy of its motif**, which is what a repeat-tract
                // genotype is about: the record has to carry the motif and both flanks, and the
                // call has to say how many copies each allele holds.
                let mut cigar = vec![Op::new(Kind::Match, sequence.len())];
                if carries_the_alternative
                    && shortens_the_tract
                    && start < tract_from
                    && start + sequence.len() > tract_to
                {
                    let cut = tract_to - TRACT_MOTIF.len() - start;
                    sequence.drain(cut..cut + TRACT_MOTIF.len());
                    cigar = vec![
                        Op::new(Kind::Match, cut),
                        Op::new(Kind::Deletion, TRACT_MOTIF.len()),
                        Op::new(Kind::Match, sequence.len() - cut),
                    ];
                }
                // **A wrong base in some reads**, where asked for: at an offset hashed from the
                // sample's name and the read, so that few positions collect the same wrong base
                // from more than one read — which the caller would rightly call a variant — and
                // never on the base the read's own sample varies at, which would turn an
                // alternative read back into a reference one. (The comparison is in reference
                // coordinates, which a read the tract's deletion shortened shifts past the cut;
                // no such read covers a sample's site here.) The two constants are splitmix64's
                // golden-ratio and mixing multipliers.
                if with_sequencing_errors {
                    let of_the_name = sample.bytes().fold(0_u64, |hash, byte| {
                        (hash ^ u64::from(byte)).wrapping_mul(0xD1B5_4A32_D192_ED03)
                    });
                    let mut hash =
                        (index as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ of_the_name;
                    hash ^= hash >> 29;
                    let offset = (hash % sequence.len() as u64) as usize;
                    // About one read in four, each to one of the three other bases, so two reads
                    // rarely err at one position and two that do rarely agree.
                    let errs = (hash >> 40).is_multiple_of(4);
                    if errs && start + offset != substitution {
                        let mut others = [0_u8; 3];
                        for (slot, base) in b"ACGT"
                            .iter()
                            .filter(|base| **base != sequence[offset])
                            .enumerate()
                        {
                            others[slot] = *base;
                        }
                        sequence[offset] = others[((hash >> 20) % 3) as usize];
                    }
                }
                // **Both strands, and the alternative reads lean to one of them.** The
                // artifact test compares the alternative reads' forward share against the
                // reference reads' own, so a fixture that split both evenly makes the two
                // shares equal — and equal is also what a stored forward count destroyed on
                // write gives, since zero over zero agrees with zero over zero. Measured: with
                // the strands balanced, zeroing the count on write left this oracle green.
                // Three alternative reads in four are forward here and one reference read in
                // four is, which is a bias the test scores and a destroyed count erases.
                let forward = if carries_the_alternative {
                    index % 8 != 6
                } else {
                    index % 8 == 1
                };
                let flags = if forward {
                    Flags::empty()
                } else {
                    Flags::REVERSE_COMPLEMENTED
                };
                let mut record = RecordBuf::builder()
                    .set_name(format!("{sample}-r{index}").as_bytes())
                    .set_reference_sequence_id(0_usize)
                    .set_flags(flags)
                    .set_mapping_quality(MappingQuality::new(60).expect("mapq in range"))
                    .set_alignment_start(
                        noodles_core::Position::try_from(start + 1).expect("a 1-based position"),
                    )
                    .set_cigar(Cigar::from(cigar))
                    .set_quality_scores(QualityScores::from(vec![35_u8; sequence.len()]))
                    .set_sequence(Sequence::from(sequence))
                    .build();
                // **Read groups alternate**, so a sample declaring two has reads in both and a
                // walk numbers two groups from zero.
                record.data_mut().insert(
                    Tag::READ_GROUP,
                    Value::String(
                        read_groups[index % read_groups.len()]
                            .as_bytes()
                            .to_vec()
                            .into(),
                    ),
                );
                record
            })
            .collect()
    };

    let file_of = |sample: &str,
                   file: &str,
                   substitution: usize,
                   shortens_the_tract: bool,
                   read_groups: &[String]| {
        let declared: Vec<FixtureReadGroup<'_>> = read_groups
            .iter()
            .map(|id| FixtureReadGroup::new(id, Some(sample)))
            .collect();
        indexed_named_bam(
            &header_with_read_groups(Some("coordinate"), &[(name, length, None)], &declared),
            &reads_of(sample, substitution, shortens_the_tract, read_groups),
            file,
        )
    };
    // **The first sample declares two read groups and the second one**, and no id is shared:
    // a run refuses a cohort in which two read groups anywhere carry one `@RG ID`.
    let first_groups = ["rg-one-a".to_owned(), "rg-one-b".to_owned()];
    let second_groups = ["rg-two".to_owned()];
    let (first, first_bam) = file_of(
        "one",
        "one.bam",
        FIRST_SAMPLES_SUBSTITUTION,
        true,
        &first_groups,
    );
    let (second, second_bam) = file_of(
        "two",
        "two.bam",
        SECOND_SAMPLES_SUBSTITUTION,
        false,
        &second_groups,
    );

    let mut samples = vec!["one".to_owned(), "two".to_owned()];
    let mut alignments = vec![first_bam, second_bam];
    let mut files = vec![first, second];
    // **Further samples' sites**: every 23rd base from the 31st, skipping the tract with a margin
    // and the first two samples' sites, so no two samples share a site.
    let mut sites = (0..)
        .map(|n| 31 + 23 * n)
        .filter(|site| {
            !(TRACT.0 - 5..TRACT.1 + 5).contains(site)
                && ![FIRST_SAMPLES_SUBSTITUTION, SECOND_SAMPLES_SUBSTITUTION].contains(site)
        })
        .take_while(|site| *site < length - 30);
    for further in 0..further_samples {
        let sample = format!("s{:02}", further + 3);
        let site = sites
            .next()
            .expect("the contig holds a site for every further sample asked for");
        let (file, bam) = file_of(
            &sample,
            &format!("{sample}.bam"),
            site,
            false,
            &[format!("rg-{sample}")],
        );
        samples.push(sample);
        alignments.push(bam);
        files.push(file);
    }

    AVaryingCohort {
        reference,
        catalog,
        alignments,
        samples,
        directory,
        _files: files,
    }
}
