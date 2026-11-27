//! **`pvc call-from-psps --regions` reaching outside the psps' ground, as a person meets it**: the
//! `pvc` binary itself, run on files it wrote, so that the exit status and the text on the terminal
//! are what is checked — not a library function's `Err`.
//!
//! The library's own tests (`src/cli/call_from_psps/tests.rs`) pin the message word for word and
//! that the refusal comes before the reference is read. This one adds what only the process can
//! show: the run ends with a failing exit status, prints the refusal, and leaves no VCF behind.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use noodles_bam as bam;
use noodles_sam as sam;

/// Run `pvc` with `arguments`, from `directory`.
fn pvc(directory: &Path, arguments: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_pvc"))
        .args(arguments)
        .current_dir(directory)
        .output()
        .expect("pvc starts")
}

/// A reference of two contigs, `chr1` of 300 bases and `chr2` of 200, of mixed bases with no
/// repeat to speak of.
fn a_reference(directory: &Path) -> PathBuf {
    let mut state: u32 = 12_345;
    let mut bases = |count: usize| -> String {
        (0..count)
            .map(|_| {
                state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                b"ACGT"[(state >> 16) as usize % 4] as char
            })
            .collect()
    };
    let path = directory.join("ref.fa");
    let mut fasta = File::create(&path).expect("the FASTA is created");
    for (name, length) in [("chr1", 300), ("chr2", 200)] {
        writeln!(fasta, ">{name}\n{}", bases(length)).expect("the FASTA is written");
    }
    path
}

/// A BAM of one sample over that reference, holding no reads: enough for a walk to write a psp.
fn an_empty_alignment_file(directory: &Path) -> PathBuf {
    let header: sam::Header = "@HD\tVN:1.6\tSO:coordinate\n\
                               @SQ\tSN:chr1\tLN:300\n\
                               @SQ\tSN:chr2\tLN:200\n\
                               @RG\tID:rg-one\tSM:one\n"
        .parse()
        .expect("the SAM header parses");
    let path = directory.join("one.bam");
    let mut writer = bam::io::Writer::new(File::create(&path).expect("the BAM is created"));
    writer.write_header(&header).expect("the header is written");
    writer.try_finish().expect("the BAM is finished");
    path
}

/// **The run stops with a failing exit status, prints which pieces are outside, and writes no
/// VCF**, when `--regions` asks for ground the psps were not written over. The psps are walked
/// over `chr1:1-100`; the run asks for `chr1:51-150`.
#[test]
fn pvc_refuses_regions_outside_the_psps_with_a_failing_exit_status() {
    let directory = tempfile::tempdir().expect("a temporary directory");
    let here = directory.path();
    let reference = a_reference(here);
    let alignments = an_empty_alignment_file(here);
    std::fs::write(here.join("walked.bed"), "chr1\t0\t100\n").expect("the walked BED");
    std::fs::write(here.join("asked.bed"), "chr1\t50\t150\n").expect("the asked BED");
    let reference = reference.to_str().expect("a UTF-8 path");
    let alignments = alignments.to_str().expect("a UTF-8 path");

    let catalog = pvc(here, &["repeat-catalog", "--reference", reference]);
    assert!(
        catalog.status.success(),
        "{}",
        String::from_utf8_lossy(&catalog.stderr)
    );
    let walk = pvc(
        here,
        &[
            "generate-psps",
            "--reference",
            reference,
            "--alignment",
            alignments,
            "--build-index-if-missing",
            "--regions",
            "walked.bed",
            "--output-dir",
            "psps",
        ],
    );
    assert!(
        walk.status.success(),
        "{}",
        String::from_utf8_lossy(&walk.stderr)
    );

    let call = pvc(
        here,
        &[
            "call-from-psps",
            "--reference",
            reference,
            "--psp",
            "psps",
            "--defaults",
            "--regions",
            "asked.bed",
            "--output",
            "calls.vcf",
        ],
    );

    let printed = String::from_utf8_lossy(&call.stderr);
    assert_eq!(call.status.code(), Some(1), "{printed}");
    for line in [
        "--regions asked.bed asks for regions the psps do not cover, so nothing was called:",
        "  chr1:101-150   on chr1 the psps cover chr1:1-100",
        "The psps cover chr1:1-100 (100 bases). Positions are 1-based, both ends \
         included.",
        "Remove these regions from the BED, or write the psps again over regions that include them.",
    ] {
        assert!(
            printed.contains(line),
            "the refusal prints {line:?}, and printed:\n{printed}"
        );
    }
    assert!(!here.join("calls.vcf").exists(), "no VCF is written");
    assert!(!here.join("calls.vcf.tmp").exists(), "nor begun");
}
