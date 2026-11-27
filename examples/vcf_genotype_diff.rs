//! **How two VCFs of the same records differ, sample by sample** — a measurement tool for changes
//! to refinement that are not meant to move the output much.
//!
//! ```text
//! vcf_genotype_diff <a.vcf> <b.vcf>
//! ```
//!
//! Both plain text, the same records in the same order (header lines skipped). Counts, over every
//! sample at every record: genotypes that differ (alleles compared unordered, so phase alone is not
//! a difference), calls made in one and missing in the other, phased in one only, and the mean
//! absolute change of `GQ` and of the largest `GP`.
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process::ExitCode;

#[derive(Default)]
struct Tally {
    records: u64,
    genotypes: u64,
    differ: u64,
    differ_hom_het: u64,
    missing_in_one: u64,
    phased_in_one: u64,
    gq_change: f64,
    gq_signed: f64,
    gq_dropped_10: u64,
    gq_counted: u64,
    gp_change: f64,
    gp_counted: u64,
}

fn alleles(gt: &str) -> Option<(u32, u32)> {
    let mut parts = gt.split(['/', '|']);
    let a = parts.next()?.parse().ok()?;
    let b = parts.next()?.parse().ok()?;
    Some(if a <= b { (a, b) } else { (b, a) })
}

fn field<'a>(keys: &[&str], values: &[&'a str], key: &str) -> Option<&'a str> {
    keys.iter()
        .position(|k| *k == key)
        .and_then(|at| values.get(at).copied())
}

fn largest(gp: &str) -> Option<f64> {
    gp.split(',')
        .map(|value| value.parse::<f64>().ok())
        .try_fold(f64::NEG_INFINITY, |best, value| value.map(|v| best.max(v)))
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("usage: vcf_genotype_diff <a.vcf> <b.vcf>");
        return ExitCode::from(2);
    }
    let open = |path: &str| BufReader::with_capacity(16 << 20, File::open(path).expect("opens"));
    let (mut a, mut b) = (open(&args[1]), open(&args[2]));
    let (mut line_a, mut line_b) = (String::new(), String::new());
    let mut tally = Tally::default();
    let (mut only_a, mut only_b) = (0_u64, 0_u64);
    // A record either file dropped (no alternative allele left) is skipped: the other file is
    // advanced past it by position, then by its REF and ALT.
    let key = |line: &str| -> (u64, String) {
        let mut columns = line.split('\t');
        let _ = columns.next();
        let position = columns.next().and_then(|p| p.parse().ok()).unwrap_or(0);
        let _ = columns.next();
        let rest = format!(
            "{}\t{}",
            columns.next().unwrap_or(""),
            columns.next().unwrap_or("")
        );
        (position, rest)
    };
    let next = |reader: &mut BufReader<File>, line: &mut String| -> bool {
        loop {
            line.clear();
            if reader.read_line(line).expect("reads") == 0 {
                return false;
            }
            if !line.starts_with('#') {
                return true;
            }
        }
    };
    let (mut more_a, mut more_b) = (next(&mut a, &mut line_a), next(&mut b, &mut line_b));
    while more_a && more_b {
        let (key_a, key_b) = (key(&line_a), key(&line_b));
        if key_a != key_b {
            if key_a.0 < key_b.0 || (key_a.0 == key_b.0 && key_a.1 < key_b.1) {
                only_a += 1;
                more_a = next(&mut a, &mut line_a);
            } else {
                only_b += 1;
                more_b = next(&mut b, &mut line_b);
            }
            continue;
        }
        let columns_a: Vec<&str> = line_a.trim_end().split('\t').collect();
        let columns_b: Vec<&str> = line_b.trim_end().split('\t').collect();
        tally.records += 1;
        let keys_a: Vec<&str> = columns_a[8].split(':').collect();
        let keys_b: Vec<&str> = columns_b[8].split(':').collect();
        for (sample_a, sample_b) in columns_a[9..].iter().zip(&columns_b[9..]) {
            let values_a: Vec<&str> = sample_a.split(':').collect();
            let values_b: Vec<&str> = sample_b.split(':').collect();
            tally.genotypes += 1;
            let gt_a = field(&keys_a, &values_a, "GT").unwrap_or(".");
            let gt_b = field(&keys_b, &values_b, "GT").unwrap_or(".");
            match (alleles(gt_a), alleles(gt_b)) {
                (Some(x), Some(y)) => {
                    if x != y {
                        tally.differ += 1;
                        if (x.0 == x.1) != (y.0 == y.1) {
                            tally.differ_hom_het += 1;
                        }
                    }
                }
                (None, None) => {}
                _ => tally.missing_in_one += 1,
            }
            if gt_a.contains('|') != gt_b.contains('|') {
                tally.phased_in_one += 1;
            }
            if let (Some(q_a), Some(q_b)) = (
                field(&keys_a, &values_a, "GQ").and_then(|q| q.parse::<f64>().ok()),
                field(&keys_b, &values_b, "GQ").and_then(|q| q.parse::<f64>().ok()),
            ) {
                tally.gq_change += (q_a - q_b).abs();
                tally.gq_signed += q_b - q_a;
                if q_a - q_b >= 10.0 {
                    tally.gq_dropped_10 += 1;
                }
                tally.gq_counted += 1;
            }
            if let (Some(p_a), Some(p_b)) = (
                field(&keys_a, &values_a, "GP").and_then(largest),
                field(&keys_b, &values_b, "GP").and_then(largest),
            ) {
                tally.gp_change += (p_a - p_b).abs();
                tally.gp_counted += 1;
            }
        }
        more_a = next(&mut a, &mut line_a);
        more_b = next(&mut b, &mut line_b);
    }
    while more_a {
        only_a += 1;
        more_a = next(&mut a, &mut line_a);
    }
    while more_b {
        only_b += 1;
        more_b = next(&mut b, &mut line_b);
    }
    println!("records only in the first {only_a}, only in the second {only_b}");
    println!(
        "mean GQ second less first {:.4}; GQ down by 10 or more {} ({:.4}%)",
        tally.gq_signed / tally.gq_counted.max(1) as f64,
        tally.gq_dropped_10,
        100.0 * tally.gq_dropped_10 as f64 / tally.gq_counted.max(1) as f64
    );
    let share = |count: u64| 100.0 * count as f64 / tally.genotypes.max(1) as f64;
    println!(
        "records {}, genotypes {}: differ {} ({:.4}%), of which hom/het {} ({:.4}%); missing in one {} ({:.4}%); phased in one {} ({:.4}%); mean |dGQ| {:.4}; mean |d max GP| {:.6}",
        tally.records,
        tally.genotypes,
        tally.differ,
        share(tally.differ),
        tally.differ_hom_het,
        share(tally.differ_hom_het),
        tally.missing_in_one,
        share(tally.missing_in_one),
        tally.phased_in_one,
        share(tally.phased_in_one),
        tally.gq_change / tally.gq_counted.max(1) as f64,
        tally.gp_change / tally.gp_counted.max(1) as f64,
    );
    ExitCode::SUCCESS
}
