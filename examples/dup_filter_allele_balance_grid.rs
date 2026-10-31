//! **Where does the hidden-duplication filter's score put an ordinary heterozygote, as a function
//! of read depth?**
//!
//! Calls the real scorer, [`score_locus_for_paralogy`], on one sample at a time over a grid of
//! read depth, ALT share and relative copy number, with the inbreeding coefficient at zero (the
//! GIAB samples are outbred, and `--defaults` carries `F = 0`). For each cell it prints the
//! log-likelihood ratio, the posterior a run with prior `π` would give, and whether a run whose
//! cut is `LR ≥ lr_cut` would remove it.
//!
//! `π`, the cut and σ₀ are a real run's: read them from the `##paralogFilter=` header line and the
//! run log's coverage-model line, and pass them in. The defaults are placeholders.
//!
//! ```text
//! cargo run --release --example dup_filter_allele_balance_grid -- \
//!     --sigma0 0.14 --pi 0.03 --lr-cut 4.2 [--dispersion 0.0015] > grid.tsv
//! cargo run --release --example dup_filter_allele_balance_grid -- \
//!     --sigma0 0.14 --pi 0.03 --lr-cut 4.2 --site 59,245,1.10 --site 58,314,1.2
//! ```
//!
//! `--site ALT,TOTAL,COPY_NUMBER` scores single observations instead of the grid. `--dispersion` sets the
//! read-split dispersion ρ under both stories (default the run's, 0.0015; zero is production's
//! strict binomial).

use pop_var_caller::float;
use pop_var_caller::paralog::{
    DEFAULT_ALLELE_DISPERSION, LocusObservations, ParalogModelParams, ParalogScorePrecompute,
    SampleObservation, score_locus_for_paralogy,
};

struct Options {
    sigma0: f64,
    pi: f64,
    lr_cut: f64,
    dispersion: f64,
    sites: Vec<(u32, u32, f64)>,
}

fn parse_options() -> Options {
    let mut options = Options {
        sigma0: 0.14,
        pi: 0.03,
        lr_cut: 4.2,
        dispersion: DEFAULT_ALLELE_DISPERSION,
        sites: Vec::new(),
    };
    let mut args = std::env::args().skip(1);
    while let Some(flag) = args.next() {
        let value = args
            .next()
            .unwrap_or_else(|| panic!("{flag} needs a value"));
        match flag.as_str() {
            "--sigma0" => options.sigma0 = value.parse().expect("--sigma0 is a number"),
            "--pi" => options.pi = value.parse().expect("--pi is a number"),
            "--lr-cut" => options.lr_cut = value.parse().expect("--lr-cut is a number"),
            "--dispersion" => {
                options.dispersion = value.parse().expect("--dispersion is a number");
            }
            "--site" => {
                let fields: Vec<&str> = value.split(',').collect();
                assert_eq!(fields.len(), 3, "--site is ALT,TOTAL,COPY_NUMBER");
                options.sites.push((
                    fields[0].parse().expect("ALT is a count"),
                    fields[1].parse().expect("TOTAL is a count"),
                    fields[2].parse().expect("COPY_NUMBER is a number"),
                ));
            }
            other => panic!("unknown option {other}"),
        }
    }
    options
}

/// One sample's ratio, through the real scorer with production's default model constants.
fn ratio(
    precompute: &ParalogScorePrecompute,
    sigma0: f64,
    alt_reads: u32,
    total_reads: u32,
    relative_copy_number: f64,
) -> f64 {
    let samples = [Some(SampleObservation {
        relative_copy_number,
        alt_reads,
        total_reads,
        inbreeding_coefficient: 0.0,
    })];
    score_locus_for_paralogy(
        &LocusObservations { samples: &samples },
        &[sigma0],
        precompute,
    )
    .paralog_log_likelihood_ratio
}

fn posterior(lr: f64, pi: f64) -> f64 {
    1.0 / (1.0 + float::exp(-(lr + float::ln(pi / (1.0 - pi)))))
}

fn main() {
    let options = parse_options();
    let params = ParalogModelParams {
        allele_dispersion: options.dispersion,
        ..ParalogModelParams::default()
    };
    let precompute = ParalogScorePrecompute::new(&params, &[0.0]);

    println!("alt_reads\ttotal_reads\talt_share\trelative_copy_number\tlr\tposterior\tremoved");
    let emit = |alt: u32, total: u32, copy_number: f64| {
        let lr = ratio(&precompute, options.sigma0, alt, total, copy_number);
        println!(
            "{alt}\t{total}\t{:.3}\t{copy_number:.2}\t{lr:.3}\t{:.4}\t{}",
            f64::from(alt) / f64::from(total),
            posterior(lr, options.pi),
            u8::from(lr >= options.lr_cut)
        );
    };

    if !options.sites.is_empty() {
        for &(alt, total, copy_number) in &options.sites {
            emit(alt, total, copy_number);
        }
        return;
    }

    const DEPTHS: [u32; 14] = [3, 5, 8, 10, 15, 20, 30, 50, 75, 100, 150, 200, 300, 400];
    const COPY_NUMBERS: [f64; 8] = [0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5];
    for &total in &DEPTHS {
        for share_step in 0..=12 {
            // ALT share 0.20, 0.25, …, 0.80.
            let share = 0.20 + 0.05 * f64::from(share_step);
            let alt = (share * f64::from(total)).round() as u32;
            for &copy_number in &COPY_NUMBERS {
                emit(alt, total, copy_number);
            }
        }
    }
}
