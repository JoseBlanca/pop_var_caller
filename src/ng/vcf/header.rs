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

/// The VCF version every ng file declares.
///
/// **4.4**, matching both production writers. The version is what tells a reader which spelling
/// rules apply — among them the padding rule for an allele at a contig's first base, which spec
/// §5 relies on.
pub const FILE_FORMAT: &str = "VCFv4.4";

/// The nine fixed column headings, before the sample names.
const FIXED_COLUMN_HEADINGS: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

/// **Every `INFO` field the file can carry, declared** — spec §6, in the order records write
/// them.
const INFO_DECLARATIONS: &[&str] = &[
    r#"##INFO=<ID=AF,Number=A,Type=Float,Description="Fitted frequency of each ALT allele, from the calling loop's converged pass">"#,
    r#"##INFO=<ID=AC,Number=A,Type=Integer,Description="Copies of each ALT allele in the called genotypes">"#,
    r#"##INFO=<ID=AN,Number=1,Type=Integer,Description="Total called allele copies (no-call samples excluded)">"#,
    r#"##INFO=<ID=DP,Number=1,Type=Integer,Description="Sum of the samples' DP">"#,
    r#"##INFO=<ID=ABPEN,Number=1,Type=Float,Description="Phred subtracted from QUAL by the allele-balance artifact test">"#,
    r#"##INFO=<ID=SPPEN,Number=1,Type=Float,Description="Phred subtracted from QUAL by the strand and read-position artifact test">"#,
    r#"##INFO=<ID=MQREF,Number=1,Type=Float,Description="Cohort-pooled mean mapping quality of reads supporting REF">"#,
    r#"##INFO=<ID=MQALT,Number=A,Type=Float,Description="Cohort-pooled mean mapping quality of reads supporting each ALT">"#,
    r#"##INFO=<ID=MQDIFF,Number=A,Type=Float,Description="MQALT minus MQREF per ALT; negative means ALT reads map worse (multi-mapper fingerprint)">"#,
    r#"##INFO=<ID=STR,Number=0,Type=Flag,Description="This record is a repeat-tract locus">"#,
    r#"##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit of the tract, reference strand">"#,
    r#"##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Repeat unit length in bases">"#,
];

/// **Every `FILTER` value the file can carry, declared** — spec §8, `PASS` included.
///
/// Production's repeat-tract writer leaves `PASS` undeclared, which is legal VCF and gratuitous:
/// declaring it costs one line and spares a reader wondering whether the file means something
/// else by it.
const FILTER_DECLARATIONS: &[&str] = &[
    r#"##FILTER=<ID=PASS,Description="All filters passed">"#,
    r#"##FILTER=<ID=EMNoConv,Description="The calling loop did not converge within its pass cap">"#,
    r#"##FILTER=<ID=notPeriodic,Description="Tract allele-length distribution inconsistent with the motif period">"#,
    r#"##FILTER=<ID=tooManyAlleles,Description="More candidate alleles segregate than the caller admits">"#,
    r#"##FILTER=<ID=lowDepth,Description="Insufficient cohort depth to call the tract">"#,
];

/// **Every `FORMAT` field the file can carry, declared** — spec §7.
const FORMAT_DECLARATIONS: &[&str] = &[
    r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype, unphased">"#,
    r#"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Phred probability the called genotype is wrong, capped at 99">"#,
    r#"##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Reads this sample observed at the locus, whether or not a written allele explains them">"#,
    r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Reads whose sequence matched each allele exactly, REF first">"#,
    r#"##FORMAT=<ID=REPCN,Number=.,Type=Integer,Description="Repeat copy number of each called allele, GT order">"#,
];

/// **The whole header, as the file's text**, ending with a newline so the first record follows
/// directly.
///
/// The order is spec §4's: the format version, what wrote the file and how, what it was called
/// against, the contigs, then the declarations, then the column headings.
///
/// **Hand-written rather than assembled through noodles**, which the plan proposed. noodles
/// groups header records by kind in an order of its own, and §4 fixes a different one; since the
/// records are hand-written text already, one mechanism and one order is worth more here than
/// the library's validation, and `bcftools` is what judges the result.
///
/// **A provenance line whose value is empty is omitted rather than written empty**, which is
/// production's rule: `##reference=` states nothing and invites a reader to think the run had a
/// reference it could not name.
#[must_use]
pub fn header_text(metadata: &VcfHeaderMetadata) -> String {
    let mut lines: Vec<String> = vec![format!("##fileformat={FILE_FORMAT}")];

    for (key, value) in [
        ("source", metadata.source()),
        ("commandline", metadata.command_line().to_string()),
        ("reference", metadata.reference_path().to_string()),
        (
            "parametersFile",
            metadata.parameters_file_name().to_string(),
        ),
    ] {
        if !value.is_empty() {
            lines.push(format!("##{key}={value}"));
        }
    }

    for contig in metadata.contigs() {
        lines.push(contig_line(contig));
    }

    lines.extend(INFO_DECLARATIONS.iter().map(ToString::to_string));
    lines.extend(FILTER_DECLARATIONS.iter().map(ToString::to_string));
    lines.extend(FORMAT_DECLARATIONS.iter().map(ToString::to_string));

    let mut headings = FIXED_COLUMN_HEADINGS.to_string();
    for name in metadata.sample_names() {
        headings.push('\t');
        headings.push_str(name);
    }
    lines.push(headings);

    let mut text = lines.join("\n");
    text.push('\n');
    text
}

/// One `##contig` line. The digest is written only where the run has one — see
/// [`HeaderContig::md5`].
fn contig_line(contig: &HeaderContig) -> String {
    let mut line = format!("##contig=<ID={},length={}", contig.name, contig.length);
    if let Some(digest) = contig.md5 {
        line.push_str(",md5=");
        for byte in digest {
            line.push_str(&format!("{byte:02x}"));
        }
    }
    line.push('>');
    line
}

#[cfg(test)]
mod header_text_tests {
    use super::*;

    fn metadata_with(contigs: Vec<HeaderContig>, samples: &[&str]) -> VcfHeaderMetadata {
        VcfHeaderMetadata::try_new(
            contigs,
            samples.iter().map(|name| (*name).to_string()).collect(),
            "ng call --reference ref.fa cohort/*.cram".to_string(),
            "/genomes/ref.fa".to_string(),
            "run.parameters.toml".to_string(),
        )
        .expect("a well-formed header")
    }

    #[test]
    fn the_whole_header_of_a_two_sample_run() {
        let metadata = metadata_with(
            vec![HeaderContig {
                name: "chr1".to_string(),
                length: 248_956_422,
                md5: None,
            }],
            &["HG002", "HG003"],
        );

        let expected = format!(
            "##fileformat=VCFv4.4\n\
             ##source=ng {version}\n\
             ##commandline=ng call --reference ref.fa cohort/*.cram\n\
             ##reference=/genomes/ref.fa\n\
             ##parametersFile=run.parameters.toml\n\
             ##contig=<ID=chr1,length=248956422>\n\
             {info}\n{filter}\n{format}\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\tHG003\n",
            version = env!("CARGO_PKG_VERSION"),
            info = INFO_DECLARATIONS.join("\n"),
            filter = FILTER_DECLARATIONS.join("\n"),
            format = FORMAT_DECLARATIONS.join("\n"),
        );

        assert_eq!(header_text(&metadata), expected);
    }

    #[test]
    fn a_contig_digest_is_written_as_hexadecimal_when_the_run_read_the_bases() {
        let metadata = metadata_with(
            vec![HeaderContig {
                name: "chr1".to_string(),
                length: 1_000,
                md5: Some([
                    0x6a, 0xef, 0x89, 0x7c, 0x3d, 0x6f, 0xf0, 0xc7, 0x8a, 0xff, 0x06, 0xac, 0x18,
                    0x91, 0x78, 0xdd,
                ]),
            }],
            &["one"],
        );

        assert!(
            header_text(&metadata)
                .contains("##contig=<ID=chr1,length=1000,md5=6aef897c3d6ff0c78aff06ac189178dd>"),
            "got {}",
            header_text(&metadata)
        );
    }

    #[test]
    fn a_run_without_the_reference_bases_writes_a_contig_line_with_no_digest() {
        let metadata = metadata_with(
            vec![HeaderContig {
                name: "chr1".to_string(),
                length: 1_000,
                md5: None,
            }],
            &["one"],
        );
        assert!(header_text(&metadata).contains("##contig=<ID=chr1,length=1000>\n"));
        assert!(!header_text(&metadata).contains("md5="));
    }

    #[test]
    fn a_provenance_line_with_nothing_to_say_is_omitted_rather_than_written_empty() {
        // Production's rule: `##reference=` states nothing and invites a reader to think the
        // run had a reference it could not name.
        let metadata = VcfHeaderMetadata::try_new(
            Vec::new(),
            vec!["one".to_string()],
            String::new(),
            String::new(),
            "run.parameters.toml".to_string(),
        )
        .expect("a well-formed header");

        let text = header_text(&metadata);
        assert!(!text.contains("##commandline"), "got {text}");
        assert!(!text.contains("##reference"), "got {text}");
        // The one that is present still is.
        assert!(text.contains("##parametersFile=run.parameters.toml\n"));
        // And `##source` is derived, so it is always there.
        assert!(text.contains("##source=ng "));
    }

    #[test]
    fn every_value_the_records_can_write_is_declared() {
        // The rule that keeps a file self-describing: nothing may appear in a record that the
        // header did not declare. These are the ids the encoder emits.
        let metadata = metadata_with(Vec::new(), &["one"]);
        let text = header_text(&metadata);

        for id in [
            "AF", "AC", "AN", "DP", "ABPEN", "SPPEN", "MQREF", "MQALT", "MQDIFF", "STR", "RU",
            "PERIOD",
        ] {
            assert!(text.contains(&format!("##INFO=<ID={id},")), "INFO {id}");
        }
        for id in ["GT", "GQ", "DP", "AD", "REPCN"] {
            assert!(text.contains(&format!("##FORMAT=<ID={id},")), "FORMAT {id}");
        }
        for id in [
            "PASS",
            "EMNoConv",
            "notPeriodic",
            "tooManyAlleles",
            "lowDepth",
        ] {
            assert!(text.contains(&format!("##FILTER=<ID={id},")), "FILTER {id}");
        }
    }

    #[test]
    fn every_filter_verdict_the_encoder_can_write_has_a_declaration() {
        // Stronger than the list above: it walks the enum, so a verdict added later without a
        // declaration fails here rather than in a consumer.
        use crate::ng::vcf::FilterVerdict;

        let text = header_text(&metadata_with(Vec::new(), &["one"]));
        for verdict in [
            FilterVerdict::Pass,
            FilterVerdict::EmDidNotConverge,
            FilterVerdict::NotPeriodic,
            FilterVerdict::TooManyAlleles,
            FilterVerdict::LowDepth,
        ] {
            assert!(
                text.contains(&format!("##FILTER=<ID={},", verdict.as_str())),
                "no declaration for {}",
                verdict.as_str()
            );
        }
    }

    #[test]
    fn the_header_ends_with_a_newline_so_the_first_record_follows_directly() {
        let text = header_text(&metadata_with(Vec::new(), &["one"]));
        assert!(text.ends_with('\n'));
        assert!(!text.ends_with("\n\n"));
    }

    #[test]
    fn the_sample_names_close_the_column_headings_in_the_run_s_order() {
        let text = header_text(&metadata_with(Vec::new(), &["c", "a", "b"]));
        let headings = text
            .lines()
            .next_back()
            .expect("the column headings are the last line");
        assert!(headings.ends_with("\tFORMAT\tc\ta\tb"), "got {headings}");
    }
}
