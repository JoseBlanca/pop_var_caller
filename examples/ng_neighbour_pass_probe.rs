//! **The neighbour-evidence pass's window code run over a finished VCF** — a research probe for
//! phase two of `doc/devel/ng/spec/neighbour_evidence_pass.md`, before the pass is wired into a
//! run.
//!
//! ```text
//! ng_neighbour_pass_probe <calls.vcf> <parameters.toml> <out.tsv> [core_kb margin_kb floor_kb]
//! ```
//!
//! Reads an uncompressed pvc VCF and the run's parameters file, keeps the PASS records with rows,
//! plans each contig's windows and solves them, and writes one line per sample per record:
//! `contig pos ref alt sample row calling_gt outcome pass_gt pass_probability`. `row` is
//! `informative`, `flat` or `absent`; `outcome` is `unchanged`, `called` or `withheld`. Scoring
//! against a truth set is the driver's job. Window sizes default to the spec's; smaller ones make a
//! small simulated genome hold several windows.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process::ExitCode;

use pop_var_caller::neighbour_evidence::fit::SampleFit;
use pop_var_caller::neighbour_evidence::records::{SampleColumns, SampleRow, SiteRecord};
use pop_var_caller::neighbour_evidence::refine::{Outcome, refine_window};
use pop_var_caller::neighbour_evidence::templates::MAX_TEMPLATES;
use pop_var_caller::neighbour_evidence::windows::{WindowSizes, plan_window};

/// One PASS record with rows: what is written next to the pass's result.
struct Kept {
    contig: String,
    reference: String,
    alternatives: String,
    calls: Vec<String>,
    record: SiteRecord,
}

fn genotype_text(index: usize) -> String {
    let mut seen = 0;
    for second in 0.. {
        for first in 0..=second {
            if seen == index {
                return format!("{first}/{second}");
            }
            seen += 1;
        }
    }
    unreachable!()
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 && args.len() != 7 {
        eprintln!(
            "usage: {} <calls.vcf> <parameters.toml> <out.tsv> [core_kb margin_kb floor_kb]",
            args[0]
        );
        return ExitCode::from(2);
    }
    let mut sizes = WindowSizes::default();
    if args.len() == 7 {
        let kb = |text: &str| text.parse::<u64>().expect("a whole number of kb") * 1_000;
        sizes.core_bases = kb(&args[4]);
        sizes.margin_bases = kb(&args[5]);
        sizes.margin_floor_bases = kb(&args[6]);
    }

    let parameters: toml::Table =
        toml::from_str(&std::fs::read_to_string(&args[2]).expect("the parameters file"))
            .expect("the parameters file is TOML");
    let inbreeding: HashMap<String, f64> = parameters["inbreeding"]["by_sample"]
        .as_array()
        .expect("inbreeding.by_sample")
        .iter()
        .map(|row| {
            (
                row["sample"].as_str().expect("a sample").to_string(),
                row["inbreeding_coefficient"]["value"]
                    .as_float()
                    .expect("a coefficient"),
            )
        })
        .collect();

    let input = BufReader::new(std::fs::File::open(&args[1]).expect("the VCF"));
    let mut samples: Vec<String> = Vec::new();
    let mut contig_lengths: HashMap<String, u64> = HashMap::new();
    let mut contigs: Vec<(String, Vec<Kept>)> = Vec::new();
    for line in input.lines() {
        let line = line.expect("a line");
        if let Some(contig) = line.strip_prefix("##contig=<ID=") {
            let name = contig.split(',').next().expect("an ID").to_string();
            if let Some(length) = contig.split("length=").nth(1) {
                let length = length
                    .trim_end_matches('>')
                    .split(',')
                    .next()
                    .expect("a length");
                contig_lengths.insert(name, length.parse().expect("a contig length"));
            }
            continue;
        }
        if line.starts_with("##") {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if line.starts_with('#') {
            samples = fields[9..].iter().map(|name| (*name).to_string()).collect();
            continue;
        }
        if fields[6] != "PASS" || fields[4] == "." {
            continue;
        }
        let format: Vec<&str> = fields[8].split(':').collect();
        let (Some(pl), Some(gp)) = (
            format.iter().position(|key| *key == "PL"),
            format.iter().position(|key| *key == "GP"),
        ) else {
            continue;
        };
        let alt_frequencies: Vec<f64> = fields[7]
            .split(';')
            .find_map(|entry| entry.strip_prefix("AF="))
            .expect("an AF")
            .split(',')
            .map(|value| value.parse().expect("an allele frequency"))
            .collect();
        let mut calls = Vec::with_capacity(samples.len());
        let parsed: Vec<Option<(Vec<u32>, Vec<f64>)>> = fields[9..]
            .iter()
            .map(|column| {
                let values: Vec<&str> = column.split(':').collect();
                calls.push(values[0].to_string());
                let phred = values.get(pl).copied().unwrap_or(".");
                let posteriors = values.get(gp).copied().unwrap_or(".");
                if phred == "." || posteriors == "." {
                    return None;
                }
                Some((
                    phred
                        .split(',')
                        .map(|value| value.parse().expect("a PL"))
                        .collect(),
                    posteriors
                        .split(',')
                        .map(|value| value.parse().expect("a GP"))
                        .collect(),
                ))
            })
            .collect();
        if parsed.iter().all(Option::is_none) {
            continue;
        }
        let columns: Vec<Option<SampleColumns<'_>>> = parsed
            .iter()
            .map(|column| {
                column.as_ref().map(|(phred, posteriors)| SampleColumns {
                    phred_likelihoods: phred,
                    genotype_posteriors: posteriors,
                })
            })
            .collect();
        let record = SiteRecord::new(
            fields[1].parse().expect("a POS"),
            &alt_frequencies,
            &columns,
        );
        if contigs.last().is_none_or(|(name, _)| name != fields[0]) {
            contigs.push((fields[0].to_string(), Vec::new()));
        }
        contigs.last_mut().expect("a contig").1.push(Kept {
            contig: fields[0].to_string(),
            reference: fields[3].to_string(),
            alternatives: fields[4].to_string(),
            calls,
            record,
        });
    }

    let mut fits: Vec<SampleFit> = samples
        .iter()
        .map(|sample| SampleFit::new(*inbreeding.get(sample).unwrap_or(&0.0)))
        .collect();
    let mut out = BufWriter::new(std::fs::File::create(&args[3]).expect("the output"));
    writeln!(
        out,
        "contig\tpos\tref\talt\tsample\trow\tcalling_gt\toutcome\tpass_gt\tpass_probability"
    )
    .expect("write");
    let mut windows = 0;
    let mut over_budget = 0;
    for (contig, kept) in &contigs {
        let records: Vec<SiteRecord> = kept.iter().map(|kept| kept.record.clone()).collect();
        let positions: Vec<u64> = records.iter().map(|record| record.position).collect();
        let mut core_start = 0;
        while core_start < records.len() {
            let plan = plan_window(&positions, core_start, true, samples.len(), &sizes)
                .expect("a whole contig always plans");
            windows += 1;
            over_budget = over_budget.max(plan.over_budget);
            let refined = refine_window(
                &records[plan.context_start..plan.context_end],
                plan.core_start - plan.context_start..plan.core_end - plan.context_start,
                contig_lengths.get(contig).copied(),
                &mut fits,
                MAX_TEMPLATES,
            );
            for (offset, result) in refined.iter().enumerate() {
                let kept = &kept[plan.core_start + offset];
                let genotypes = kept.record.genotype_count();
                for (sample, name) in samples.iter().enumerate() {
                    let row = match kept.record.row(sample) {
                        SampleRow::Informative => "informative",
                        SampleRow::Flat => "flat",
                        SampleRow::Absent => "absent",
                    };
                    let outcome = result.outcomes[sample];
                    let (outcome_text, pass_gt, probability) = if outcome == Outcome::Unchanged {
                        ("unchanged", ".".to_string(), f64::NAN)
                    } else {
                        let posterior =
                            &result.posteriors[sample * genotypes..(sample + 1) * genotypes];
                        let best = (0..genotypes).fold(0, |best, genotype| {
                            if posterior[genotype] > posterior[best] {
                                genotype
                            } else {
                                best
                            }
                        });
                        let text = if outcome == Outcome::Called {
                            "called"
                        } else {
                            "withheld"
                        };
                        (text, genotype_text(best), posterior[best])
                    };
                    writeln!(
                        out,
                        "{}\t{}\t{}\t{}\t{name}\t{row}\t{}\t{outcome_text}\t{pass_gt}\t{probability:.4}",
                        kept.contig,
                        kept.record.position,
                        kept.reference,
                        kept.alternatives,
                        kept.calls[sample]
                    )
                    .expect("write");
                }
            }
            core_start = plan.core_end;
        }
    }
    let frozen = fits.iter().filter(|fit| fit.is_frozen()).count();
    eprintln!(
        "{} samples, {} contigs, {windows} windows, largest over-budget {over_budget}, {frozen} fits frozen",
        samples.len(),
        contigs.len()
    );
    for (sample, fit) in samples.iter().zip(&fits) {
        let parameters = fit.parameters(0.01);
        eprintln!(
            "{sample}\tS={}\tR={}\tF={:.3}",
            parameters.segment_bases, parameters.run_bases, parameters.autozygous_share
        );
    }
    ExitCode::SUCCESS
}
