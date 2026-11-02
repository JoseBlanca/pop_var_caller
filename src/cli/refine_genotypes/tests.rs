use super::*;
use crate::neighbour_evidence::test_cohort::{CohortShape, test_cohort};
use crate::vcf::header::{
    PHASE_SET_DECLARATION, PHASED_GENOTYPE_DECLARATION, UNPHASED_GENOTYPE_DECLARATION,
};

const SAMPLES: usize = 20;

/// A VCF of the drawn cohort as calling would write it with the pass off, its parameters file, in
/// `directory`, and how many of its records some sample truly carries the alternative at.
fn calling_output(directory: &std::path::Path) -> (PathBuf, PathBuf, usize) {
    let cohort = test_cohort(CohortShape {
        samples: SAMPLES,
        founders: 4,
        records: 1_500,
        spacing: 600,
        segment_bases: 300_000.0,
        inbred: true,
        depth: 1.0,
        read_error: 0.01,
        seed: 13,
        silent_from: None,
    });
    let names: Vec<String> = (0..SAMPLES).map(|sample| format!("s{sample}")).collect();
    let mut text = String::from(
        "##fileformat=VCFv4.4\n##contig=<ID=chr1,length=1000000>\n\
         ##neighbourEvidence=<ran=no,why=\"switched off\">\n\
         ##INFO=<ID=AF,Number=A,Type=Float,Description=\"x\">\n",
    );
    text.push_str(CALLING_GP);
    text.push_str("\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
    text.push_str(&names.join("\t"));
    text.push('\n');
    for line in &cohort.lines {
        text.push_str(line);
        text.push('\n');
    }
    let vcf = directory.join("calls.vcf");
    std::fs::write(&vcf, text).expect("the VCF is written");
    let mut parameters = String::from("ploidy = 2\n\n[inbreeding]\nby_sample = [\n");
    for name in &names {
        parameters.push_str(&format!(
            "    {{ sample = \"{name}\", inbreeding_coefficient = {{ value = 0.95 }} }},\n"
        ));
    }
    parameters.push_str("]\n");
    let toml = directory.join("run.parameters.toml");
    std::fs::write(&toml, parameters).expect("the parameters are written");
    let carried = cohort
        .truth
        .iter()
        .filter(|genotypes| genotypes.iter().any(|&genotype| genotype > 0))
        .count();
    (vcf, toml, carried)
}

fn args(vcf: PathBuf, parameters: PathBuf, output: PathBuf) -> RefineGenotypesArgs {
    RefineGenotypesArgs {
        vcf,
        parameters,
        log: Some(output.with_extension("log")),
        output,
        threads: 0,
        force: true,
        phase_min_path_agreement: crate::neighbour_evidence::phase::PathAgreement::ALL,
    }
}

/// **A VCF calling wrote is refined, and the refined one is refused**: the output's header says
/// the pass ran, in place of the input's `ran=no`, with `GP`'s description changed; the records
/// it keeps are about those some sample truly carries the alternative at — the drawn cohort writes
/// a record at every drawn site, and at most of the rare ones no founder carries the alternative,
/// so the pass drops them as calling would have; and feeding the output back is an error that
/// names why (§9 Q4).
#[test]
fn a_calling_vcf_is_refined_and_a_refined_one_refused() {
    let directory = tempfile::tempdir_in("tmp").expect("a directory to work in");
    let (vcf, parameters, carried) = calling_output(directory.path());
    let output = directory.path().join("refined.vcf");
    run_refine_genotypes(&args(vcf.clone(), parameters.clone(), output.clone()))
        .expect("the VCF is refined");
    let text = std::fs::read_to_string(&output).expect("the output reads");
    let header: Vec<&str> = text
        .lines()
        .take_while(|line| line.starts_with('#'))
        .collect();
    assert_eq!(
        header
            .iter()
            .filter(|line| line.starts_with("##neighbourEvidence="))
            .count(),
        1
    );
    assert!(
        header
            .iter()
            .any(|line| line.starts_with("##neighbourEvidence=<ran=yes"))
    );
    assert!(header.contains(&REFINED_GP));
    // The pass writes phase: `GT`'s declaration says so, and `PS` is declared.
    assert!(header.contains(&PHASED_GENOTYPE_DECLARATION));
    assert!(header.contains(&PHASE_SET_DECLARATION));
    assert!(!header.contains(&UNPHASED_GENOTYPE_DECLARATION));
    let records_in = std::fs::read_to_string(&vcf)
        .expect("the input reads")
        .lines()
        .filter(|line| !line.starts_with('#'))
        .count();
    let records_out = text.lines().filter(|line| !line.starts_with('#')).count();
    assert!(
        records_out <= records_in
            && records_out * 100 >= carried * 95
            && records_out * 100 <= carried * 105,
        "{records_in} records in, {records_out} out; {carried} carry the alternative"
    );
    let no_calls_in = std::fs::read_to_string(&vcf)
        .expect("the input reads")
        .matches("./.")
        .count();
    assert!(text.matches("./.").count() * 2 < no_calls_in);

    let again = run_refine_genotypes(&args(
        output.clone(),
        parameters,
        directory.path().join("twice.vcf"),
    ));
    assert!(matches!(
        again,
        Err(RefineGenotypesCliError::AlreadyRefined { .. })
    ));
}

/// **A parameters file without a sample of the VCF is refused, naming the sample.**
#[test]
fn a_sample_missing_from_the_parameters_is_refused() {
    let directory = tempfile::tempdir_in("tmp").expect("a directory to work in");
    let (vcf, parameters, _) = calling_output(directory.path());
    let text = std::fs::read_to_string(&parameters).expect("the parameters read");
    std::fs::write(&parameters, text.replace("\"s7\"", "\"other\"")).expect("rewritten");
    let error = run_refine_genotypes(&args(vcf, parameters, directory.path().join("out.vcf")))
        .expect_err("refused");
    assert!(error.to_string().contains("`s7`"), "{error}");
}

/// **An input already phased is refused** (phase output spec §3): a header declaring `PS`, or a
/// `GT` written with `|`.
#[test]
fn a_phased_input_is_refused() {
    let directory = tempfile::tempdir_in("tmp").expect("a directory to work in");
    let (vcf, parameters, _) = calling_output(directory.path());
    let text = std::fs::read_to_string(&vcf).expect("the input reads");

    let declared = directory.path().join("declared.vcf");
    std::fs::write(
        &declared,
        text.replacen("##INFO=", &format!("{PHASE_SET_DECLARATION}\n##INFO="), 1),
    )
    .expect("written");
    let error = run_refine_genotypes(&args(
        declared,
        parameters.clone(),
        directory.path().join("out1.vcf"),
    ))
    .expect_err("refused");
    assert!(
        matches!(error, RefineGenotypesCliError::AlreadyPhased { .. }),
        "{error}"
    );

    let piped = directory.path().join("piped.vcf");
    let header_end = text.find("\n#CHROM").expect("a column line");
    let (header, records) = text.split_at(header_end);
    std::fs::write(
        &piped,
        format!("{header}{}", records.replacen("\t0/0:", "\t0|0:", 1)),
    )
    .expect("written");
    let error = run_refine_genotypes(&args(piped, parameters, directory.path().join("out2.vcf")))
        .expect_err("refused");
    assert!(
        matches!(error, RefineGenotypesCliError::AlreadyPhased { .. }),
        "{error}"
    );
}
