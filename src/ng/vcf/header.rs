//! **What the file says about the run that produced it** — the metadata behind every `##` line
//! and the `#CHROM` line (`doc/devel/ng/spec/vcf_output.md` §4).
//!
//! This is the header's *content*, not its text: rendering belongs to a later step. What lives
//! here is the set of facts a header states, and the refusals that stop a header from stating
//! something a reader could not act on.
//!
//! **Refused rather than asserted, and the distinction is deliberate.** The record type in this
//! module's parent panics on a bad record, because every one of its invariants binds two things
//! the same worker built moments apart — a violation is a wiring defect, not an input. Header
//! metadata is the opposite: sample names come from the alignment files, contigs from the
//! reference, the command line from whoever typed it. Two files naming the same sample is a
//! *run* someone can fix, so it is a `Result`, which is what production's header builder does
//! for the same reason ([`src/vcf/header.rs`](../../../vcf/header.rs)).

use std::collections::HashSet;

use thiserror::Error;

/// The largest contig length a VCF `##contig` line can carry.
///
/// **The VCF integer type is 32-bit signed**, so a longer contig cannot be written honestly;
/// production refuses it rather than truncating, and so does this. No real assembly is near it
/// — the largest human chromosome is about 249 million bases, one part in nine of the ceiling —
/// so this catches a corrupt reference index rather than a large genome.
pub const MAX_CONTIG_LENGTH: u64 = i32::MAX as u64;

/// **What the header states about the run**, before any of it is rendered as text.
///
/// Built once per run and refused if it says anything a reader could not act on. The
/// constructor is the only way in, so a metadata value in hand has already passed every check
/// in [`HeaderMetadataError`].
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct VcfHeaderMetadata {
    contigs: Vec<HeaderContig>,
    sample_names: Vec<String>,
    command_line: String,
    reference_path: String,
    parameters_file_name: String,
}

impl VcfHeaderMetadata {
    /// Gather what the header will state, checking what a header cannot honestly say.
    ///
    /// `sample_names` must be **in the run's sample order** — the same order every record's
    /// sample columns are in. Nothing here can check that (a permutation is still a list of
    /// distinct names); it is the run's to get right, and it is why the order is named in the
    /// parameter's own documentation rather than left implied.
    ///
    /// # Errors
    ///
    /// See [`HeaderMetadataError`]. An **empty contig list is accepted**, matching production:
    /// a run given a reference with no contigs has nothing to say about them, which is a
    /// strange run rather than an unwritable header.
    pub fn try_new(
        contigs: Vec<HeaderContig>,
        sample_names: Vec<String>,
        command_line: String,
        reference_path: String,
        parameters_file_name: String,
    ) -> Result<Self, HeaderMetadataError> {
        if sample_names.is_empty() {
            return Err(HeaderMetadataError::NoSamples);
        }
        let mut seen_samples = HashSet::with_capacity(sample_names.len());
        for name in &sample_names {
            if name.is_empty() {
                return Err(HeaderMetadataError::EmptySampleName);
            }
            if !seen_samples.insert(name.as_str()) {
                return Err(HeaderMetadataError::DuplicateSampleName(name.clone()));
            }
        }

        let mut seen_contigs = HashSet::with_capacity(contigs.len());
        for contig in &contigs {
            if contig.name.is_empty() {
                return Err(HeaderMetadataError::EmptyContigName);
            }
            if !seen_contigs.insert(contig.name.as_str()) {
                return Err(HeaderMetadataError::DuplicateContigName(
                    contig.name.clone(),
                ));
            }
            if contig.length > MAX_CONTIG_LENGTH {
                return Err(HeaderMetadataError::ContigTooLong {
                    name: contig.name.clone(),
                    length: contig.length,
                });
            }
        }

        Ok(Self {
            contigs,
            sample_names,
            command_line,
            reference_path,
            parameters_file_name,
        })
    }

    /// The contigs, in the reference's own order — one `##contig` line each.
    #[inline]
    #[must_use]
    pub fn contigs(&self) -> &[HeaderContig] {
        &self.contigs
    }

    /// The sample names, in the run's sample order — the tail of the `#CHROM` line, and the
    /// order every record's sample columns are in.
    #[inline]
    #[must_use]
    pub fn sample_names(&self) -> &[String] {
        &self.sample_names
    }

    /// The `##source` value: what wrote the file, and which build of it.
    ///
    /// **Derived rather than stored**, so a run cannot claim to have been written by something
    /// other than the binary that wrote it.
    #[inline]
    #[must_use]
    pub fn source(&self) -> String {
        format!("ng {}", env!("CARGO_PKG_VERSION"))
    }

    /// The `##commandline` value: the invocation, as it was typed.
    #[inline]
    #[must_use]
    pub fn command_line(&self) -> &str {
        &self.command_line
    }

    /// The `##reference` value: the reference this run was called against.
    #[inline]
    #[must_use]
    pub fn reference_path(&self) -> &str {
        &self.reference_path
    }

    /// The `##parametersFile` value: the parameters file written beside this VCF.
    ///
    /// **A file name, not a path**, and that is the point: the two travel as a directory, and an
    /// absolute path would be stale the first time the pair moved. It is the line that makes a
    /// run reproducible from its own output directory, and neither production writer has
    /// anything like it.
    #[inline]
    #[must_use]
    pub fn parameters_file_name(&self) -> &str {
        &self.parameters_file_name
    }
}

/// **One contig, as the header states it**: its name, its length, and its digest where the run
/// has one.
///
/// A projection of the reference's own [`ContigInfo`](crate::ng::reference_info::ContigInfo),
/// keeping the three things a `##contig` line carries and dropping the file geometry, which
/// says where bases live in a FASTA and means nothing in a VCF.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct HeaderContig {
    /// The contig's name — the same string every record's `CHROM` column names.
    pub name: String,
    /// Its length in bases.
    pub length: u64,
    /// Its MD5, where the run read the reference's bases.
    ///
    /// **`None` is honest, not missing.** A run driven from a `.fai` alone never saw the
    /// sequence, so it has no digest to state; the attribute is then left off the line rather
    /// than invented. Production's SNP/indel writer states it and its repeat-tract writer never
    /// does — this carries it when the run has it, which is neither of those two behaviours.
    pub md5: Option<[u8; 16]>,
}

/// What a header cannot honestly state.
///
/// Every one of these is reachable from a run's inputs rather than from a defect in this crate,
/// which is why they are errors and not assertions.
#[derive(Debug, Error, PartialEq, Eq)]
pub enum HeaderMetadataError {
    /// A VCF names its samples in the `#CHROM` line, and a cohort has at least one.
    #[error("a cohort has at least one sample, and this header names none")]
    NoSamples,

    /// **Two samples of one run carrying one name.** Every record's columns are positional, so
    /// nothing downstream could tell the two apart — the file would be ambiguous rather than
    /// wrong, which is worse.
    #[error(
        "two samples of this run are both named `{0}`: the file's sample columns are \
         positional, so nothing reading it could tell them apart"
    )]
    DuplicateSampleName(String),

    /// A sample with no name at all — a column heading that names nothing.
    #[error("a sample of this run has an empty name, so its column would head nothing")]
    EmptySampleName,

    /// Two contigs of one reference carrying one name, which makes every `CHROM` ambiguous.
    #[error(
        "two contigs of this reference are both named `{0}`, so the CHROM column could not \
         say which one a record is on"
    )]
    DuplicateContigName(String),

    /// A contig with no name at all.
    #[error("a contig of this reference has an empty name, so no record could name it")]
    EmptyContigName,

    /// **A contig longer than a VCF can state.** The format's integers are 32-bit signed, so
    /// the length is refused rather than truncated into a plausible smaller one.
    #[error(
        "contig `{name}` is {length} bases and a VCF states a contig length as a 32-bit \
         signed integer, so anything above {MAX_CONTIG_LENGTH} cannot be written honestly"
    )]
    ContigTooLong {
        /// Which contig.
        name: String,
        /// Its stated length.
        length: u64,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    fn contig(name: &str, length: u64) -> HeaderContig {
        HeaderContig {
            name: name.to_string(),
            length,
            md5: None,
        }
    }

    fn metadata(
        contigs: Vec<HeaderContig>,
        samples: &[&str],
    ) -> Result<VcfHeaderMetadata, HeaderMetadataError> {
        VcfHeaderMetadata::try_new(
            contigs,
            samples.iter().map(|name| (*name).to_string()).collect(),
            "ng call --reference ref.fa".to_string(),
            "/genomes/ref.fa".to_string(),
            "run.parameters.toml".to_string(),
        )
    }

    #[test]
    fn a_header_states_the_run_it_came_from() {
        let header = metadata(vec![contig("chr1", 248_956_422)], &["HG002", "HG003"])
            .expect("a well-formed header");

        assert_eq!(header.sample_names(), ["HG002", "HG003"]);
        assert_eq!(header.contigs()[0].name, "chr1");
        assert_eq!(header.command_line(), "ng call --reference ref.fa");
        assert_eq!(header.reference_path(), "/genomes/ref.fa");
        assert_eq!(header.parameters_file_name(), "run.parameters.toml");
        // The source names this binary and its version, so a file cannot claim another writer.
        assert!(header.source().starts_with("ng "));
        assert!(header.source().len() > "ng ".len());
    }

    #[test]
    fn a_contig_digest_is_absent_rather_than_invented_when_the_run_never_read_the_bases() {
        let from_index_alone =
            metadata(vec![contig("chr1", 1_000)], &["one"]).expect("a well-formed header");
        assert_eq!(from_index_alone.contigs()[0].md5, None);

        let from_the_fasta = metadata(
            vec![HeaderContig {
                name: "chr1".to_string(),
                length: 1_000,
                md5: Some([7u8; 16]),
            }],
            &["one"],
        )
        .expect("a well-formed header");
        assert_eq!(from_the_fasta.contigs()[0].md5, Some([7u8; 16]));
    }

    #[test]
    fn a_reference_with_no_contigs_is_accepted() {
        // Production accepts it, and a run with nothing to say about contigs is a strange run
        // rather than an unwritable header.
        let header = metadata(Vec::new(), &["one"]).expect("a header with no contig lines");
        assert!(header.contigs().is_empty());
    }

    #[test]
    fn the_sample_order_is_kept_exactly_as_given() {
        // It is the order every record's columns are in, so the header may not tidy it.
        let header = metadata(Vec::new(), &["c", "a", "b"]).expect("a well-formed header");
        assert_eq!(header.sample_names(), ["c", "a", "b"]);
    }

    #[test]
    fn a_cohort_with_no_samples_is_refused() {
        assert_eq!(
            metadata(vec![contig("chr1", 1_000)], &[]),
            Err(HeaderMetadataError::NoSamples)
        );
    }

    #[test]
    fn two_samples_of_one_name_are_refused() {
        assert_eq!(
            metadata(Vec::new(), &["HG002", "HG003", "HG002"]),
            Err(HeaderMetadataError::DuplicateSampleName(
                "HG002".to_string()
            ))
        );
    }

    #[test]
    fn a_sample_with_no_name_is_refused() {
        assert_eq!(
            metadata(Vec::new(), &["HG002", ""]),
            Err(HeaderMetadataError::EmptySampleName)
        );
    }

    #[test]
    fn two_contigs_of_one_name_are_refused() {
        assert_eq!(
            metadata(vec![contig("chr1", 1_000), contig("chr1", 2_000)], &["one"]),
            Err(HeaderMetadataError::DuplicateContigName("chr1".to_string()))
        );
    }

    #[test]
    fn a_contig_with_no_name_is_refused() {
        assert_eq!(
            metadata(vec![contig("", 1_000)], &["one"]),
            Err(HeaderMetadataError::EmptyContigName)
        );
    }

    #[test]
    fn a_contig_longer_than_the_format_can_state_is_refused() {
        assert_eq!(
            metadata(vec![contig("huge", MAX_CONTIG_LENGTH + 1)], &["one"]),
            Err(HeaderMetadataError::ContigTooLong {
                name: "huge".to_string(),
                length: MAX_CONTIG_LENGTH + 1,
            })
        );
    }

    #[test]
    fn a_contig_at_the_format_ceiling_is_accepted() {
        // The boundary is inclusive: exactly `i32::MAX` is writable.
        let header = metadata(vec![contig("big", MAX_CONTIG_LENGTH)], &["one"])
            .expect("a contig at the ceiling");
        assert_eq!(header.contigs()[0].length, MAX_CONTIG_LENGTH);
    }

    #[test]
    fn the_largest_real_chromosome_is_far_inside_the_ceiling() {
        // Human chr1 is about 249 million bases against a ceiling of about 2.15 billion, so
        // this refusal catches a corrupt index rather than a large genome.
        let human_chr1 = 248_956_422u64;
        assert!(human_chr1 < MAX_CONTIG_LENGTH);
        assert!(MAX_CONTIG_LENGTH / human_chr1 >= 8);
    }
}
