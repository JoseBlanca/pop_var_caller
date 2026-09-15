//! **`crate::float::exp`'s table-driven algorithm against the `libm` crate's `exp` and the
//! platform's: how long a call takes, and what each returns.**
//!
//! Plan `doc/devel/implementation_plans/portable_float.md`, step D1. `measure` times the three
//! over the argument ranges the parameter fit and the calling commands pass (A1 inventory), the
//! median of `trials` passes with the fastest and slowest beside it, alternating which goes first.
//! It writes the arguments and each side's outputs to `<output-dir>` as little-endian `f64`s, so
//! `scripts/float_exp_accuracy.py` can score each side against correctly rounded values and two
//! platforms' directories can be compared bit for bit.
//!
//! ```text
//! cargo run --release --example float_exp_check -- measure <output-dir> [trials]
//! ```

#![allow(
    clippy::disallowed_methods,
    reason = "examples are research tools outside the portable-float guarantee (clippy.toml)"
)]

use std::hint::black_box;
use std::path::Path;
use std::time::Instant;

use pop_var_caller::float;

/// Arguments per range: 2^19, as in `float_libm_vs_std`.
const ARGUMENTS_PER_RANGE: usize = 1 << 19;
const DEFAULT_TRIALS: usize = 21;
const WARM_UP_PASSES: usize = 20;

type Pass = fn(&[f64], &mut [f64]);

macro_rules! side {
    (|$x:ident| $body:expr) => {{
        fn pass(arguments: &[f64], out: &mut [f64]) {
            for (x, slot) in arguments.iter().zip(out.iter_mut()) {
                let $x: f64 = black_box(*x);
                *slot = $body;
            }
        }
        pass as Pass
    }};
}

/// The three implementations, in the order their columns are printed.
fn sides() -> [(&'static str, Pass); 3] {
    [
        ("platform", side!(|x| x.exp())),
        ("libm", side!(|x| libm::exp(x))),
        ("table", side!(|x| float::exp(x))),
    ]
}

/// The ranges, as (label, first, last): uniform over [first, last].
const RANGES: [(&str, f64, f64); 5] = [
    ("log-probabilities, [-700, 0]", -700.0, 0.0),
    ("log-ratios, [-20, 20]", -20.0, 20.0),
    ("past underflow, [-745, -700]", -745.0, -700.0),
    ("near overflow, [0, 709.7]", 0.0, 709.7),
    ("near zero, [-1e-3, 1e-3]", -1e-3, 1e-3),
];

/// Evenly spaced arguments in a scrambled order (the same permutation as `float_libm_vs_std`).
fn arguments(first: f64, last: f64) -> Vec<f64> {
    (0..ARGUMENTS_PER_RANGE)
        .map(|index| {
            let bits = ARGUMENTS_PER_RANGE.trailing_zeros();
            let scrambled = ((index as u32).wrapping_mul(0x9E37_79B9)
                & (ARGUMENTS_PER_RANGE as u32 - 1))
                .reverse_bits()
                >> (u32::BITS - bits);
            let spaced = f64::from(scrambled) / ARGUMENTS_PER_RANGE as f64;
            first + (last - first) * spaced
        })
        .collect()
}

fn timed_pass(side: Pass, arguments: &[f64], out: &mut [f64]) -> f64 {
    let start = Instant::now();
    side(arguments, out);
    black_box(&out);
    start.elapsed().as_nanos() as f64 / arguments.len() as f64
}

fn write_f64s(path: &Path, values: &[f64]) {
    let bytes: Vec<u8> = values
        .iter()
        .flat_map(|value| value.to_le_bytes())
        .collect();
    std::fs::write(path, bytes).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
}

fn measure(output_dir: &Path, trials: usize) {
    std::fs::create_dir_all(output_dir).expect("create the output directory");
    // Every line printed is also written to `<output-dir>/timings.tsv`, so a figure quoted from a
    // run can be traced to the directory its outputs are in.
    let mut report = vec![
        format!(
            "platform: {} {}; {ARGUMENTS_PER_RANGE} arguments a range; {trials} timed passes a side",
            std::env::consts::OS,
            std::env::consts::ARCH
        ),
        "range\tplatform ns/call median [min, max]\tlibm ns/call\ttable ns/call\t\
         table/platform\ttable/libm\ttable differs from libm\ttable differs from platform"
            .to_owned(),
    ];
    for line in &report {
        println!("{line}");
    }
    let sides = sides();
    for (range_index, (label, first, last)) in RANGES.iter().enumerate() {
        let arguments = arguments(*first, *last);
        let mut outputs = vec![vec![0.0; arguments.len()]; sides.len()];
        for _ in 0..WARM_UP_PASSES {
            for (side, out) in sides.iter().zip(outputs.iter_mut()) {
                timed_pass(side.1, &arguments, out);
            }
        }
        let mut timings = vec![Vec::with_capacity(trials); sides.len()];
        for trial in 0..trials {
            // Rotate which side goes first, so a drift in machine speed hits all three alike.
            for offset in 0..sides.len() {
                let which = (trial + offset) % sides.len();
                timings[which].push(timed_pass(sides[which].1, &arguments, &mut outputs[which]));
            }
        }
        let summaries: Vec<(f64, f64, f64)> = timings
            .into_iter()
            .map(|mut t| {
                t.sort_by(f64::total_cmp);
                (t[t.len() / 2], t[0], t[t.len() - 1])
            })
            .collect();
        let differing = |a: usize, b: usize| {
            outputs[a]
                .iter()
                .zip(&outputs[b])
                .filter(|(x, y)| x.to_bits() != y.to_bits())
                .count()
        };
        let line = format!(
            "{label}\t{:.2} [{:.2}, {:.2}]\t{:.2} [{:.2}, {:.2}]\t{:.2} [{:.2}, {:.2}]\t{:.2}\t{:.2}\t{}\t{}",
            summaries[0].0,
            summaries[0].1,
            summaries[0].2,
            summaries[1].0,
            summaries[1].1,
            summaries[1].2,
            summaries[2].0,
            summaries[2].1,
            summaries[2].2,
            summaries[2].0 / summaries[0].0,
            summaries[2].0 / summaries[1].0,
            differing(2, 1),
            differing(2, 0),
        );
        println!("{line}");
        report.push(line);
        write_f64s(
            &output_dir.join(format!("exp.{range_index}.arguments.f64")),
            &arguments,
        );
        for ((name, _), out) in sides.iter().zip(&outputs) {
            write_f64s(
                &output_dir.join(format!("exp.{range_index}.{name}.f64")),
                out,
            );
        }
    }
    std::fs::write(output_dir.join("timings.tsv"), report.join("\n") + "\n")
        .expect("write timings.tsv");
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("measure") if args.len() >= 3 => {
            let trials = args
                .get(3)
                .map(|t| t.parse().expect("trials is a number"))
                .unwrap_or(DEFAULT_TRIALS);
            measure(Path::new(&args[2]), trials);
        }
        _ => {
            eprintln!("usage: float_exp_check measure <output-dir> [trials]");
            std::process::exit(2);
        }
    }
}
