//! **The platform maths library against the `libm` crate, function by function: how much slower,
//! and how often the bits differ.**
//!
//! Rust's `f64::ln`, `exp`, `powf` and friends call the operating system's maths library, which
//! rounds differently on macOS and Linux. The `libm` crate is one Rust implementation that would
//! give both platforms the same bits. This program measures what swapping one for the other costs,
//! over the argument ranges the caller passes (plan `doc/devel/implementation_plans/
//! portable_float.md`, step A2).
//!
//! ```text
//! cargo run --release --example float_libm_vs_std -- measure <output-dir> [trials]
//! cargo run --release --example float_libm_vs_std -- compare <dir-a> <dir-b>
//! ```
//!
//! `measure` prints, for every function and argument range, the time a call takes through std
//! and through libm — the median of `trials` passes (default 21) with the fastest and slowest pass
//! beside it — and how many arguments give different bits. It writes every output of both sides to
//! `<output-dir>` as little-endian `f64`s. `compare` reads two such directories, one written on
//! each platform, and counts the arguments whose outputs differ: for std that is how often the two
//! platforms disagree today, and for libm it must be zero.
//!
//! `powi` is the one function here that is not a library call: Rust hands it to the compiler,
//! which either expands it or calls its runtime's `__powidf2`. It is measured against the same
//! square-and-multiply loop written in Rust.

use std::hint::black_box;
use std::path::Path;
use std::time::Instant;

/// Arguments per range. 2^19 keeps each output file at 4 MB and each timed pass under a few
/// milliseconds even for the slowest function.
const ARGUMENTS_PER_RANGE: usize = 1 << 19;

/// Default number of timed passes over a range, per side.
const DEFAULT_TRIALS: usize = 21;

/// Untimed passes over a range, per side, before timing starts.
const WARM_UP_PASSES: usize = 20;

/// One implementation's pass over a range: fills the output column from the two argument columns.
/// Each is its own loop, built by [`side!`], so the function under test is called directly and
/// may be inlined as it would be in the caller — an indirect call per argument would add the same
/// few nanoseconds to both sides and pull their ratio towards one.
type Pass = fn(&[f64], &[f64], &mut [f64]);

/// Build a [`Pass`] from a two-argument expression: `side!(|x, y| x.powf(y))`.
macro_rules! side {
    (|$x:ident, $y:pat_param| $body:expr) => {{
        fn pass(first: &[f64], second: &[f64], out: &mut [f64]) {
            for ((x, y), slot) in first.iter().zip(second).zip(out.iter_mut()) {
                let $x: f64 = black_box(*x);
                let $y: f64 = black_box(*y);
                *slot = $body;
            }
        }
        pass as Pass
    }};
}

/// One function under test, as the pair of implementations compared.
struct Function {
    name: &'static str,
    std_side: Pass,
    libm_side: Pass,
}

/// One range of arguments a function is evaluated over, drawn deterministically so both
/// platforms see the same inputs.
struct ArgumentRange {
    /// What the caller passes in this range, in words — printed beside the result.
    label: &'static str,
    /// The first argument of argument number `index` out of `ARGUMENTS_PER_RANGE`.
    first: fn(f64) -> f64,
    /// The second argument (the exponent, for `powf` and `powi`); unused otherwise.
    second: fn(f64) -> f64,
}

/// `x` in `[0, 1)` mapped log-uniformly onto `[low, high]` — through libm, so both platforms
/// generate the same arguments.
fn log_uniform(x: f64, low: f64, high: f64) -> f64 {
    libm::exp(libm::log(low) + x * (libm::log(high) - libm::log(low)))
}

/// The square-and-multiply loop compiler runtimes use for `__powidf2`, written in Rust.
fn powi_loop(base: f64, exponent: i32) -> f64 {
    let mut factor = base;
    let mut remaining = exponent.unsigned_abs();
    let mut product = 1.0;
    loop {
        if remaining & 1 != 0 {
            product *= factor;
        }
        remaining >>= 1;
        if remaining == 0 {
            break;
        }
        factor *= factor;
    }
    if exponent < 0 { 1.0 / product } else { product }
}

fn unused(_: f64) -> f64 {
    0.0
}

/// Every function and the argument ranges it is measured over.
fn cases() -> Vec<(Function, Vec<ArgumentRange>)> {
    vec![
        (
            // The loop and `black_box` alone, with no function in it: its time is the part of
            // every other row's time that is not the function, so a ratio of two rows' times
            // understates the ratio of the two functions' own costs by this much.
            Function {
                name: "loop overhead",
                std_side: side!(|x, _| x),
                libm_side: side!(|x, _| x),
            },
            vec![ArgumentRange {
                label: "no function: the argument copied through",
                first: |x| x,
                second: unused,
            }],
        ),
        (
            Function {
                name: "ln",
                std_side: side!(|x, _| x.ln()),
                libm_side: side!(|x, _| libm::log(x)),
            },
            vec![
                ArgumentRange {
                    label: "probabilities and the fit's rescaled products, log-uniform in [1e-300, 1]",
                    first: |x| log_uniform(x, 1e-300, 1.0),
                    second: unused,
                },
                ArgumentRange {
                    label: "probabilities near 1: 1 - d, d log-uniform in [1e-12, 0.5]",
                    first: |x| 1.0 - log_uniform(x, 1e-12, 0.5),
                    second: unused,
                },
                ArgumentRange {
                    label: "counts and concentrations, log-uniform in [1, 1e4]",
                    first: |x| log_uniform(x, 1.0, 1e4),
                    second: unused,
                },
                ArgumentRange {
                    label: "large concentrations and odds, log-uniform in [1e4, 1e12]",
                    first: |x| log_uniform(x, 1e4, 1e12),
                    second: unused,
                },
                ArgumentRange {
                    label: "probabilities at the floor, log-uniform in [2.2e-308, 1e-300]",
                    first: |x| log_uniform(x, f64::MIN_POSITIVE, 1e-300),
                    second: unused,
                },
            ],
        ),
        (
            Function {
                name: "exp",
                std_side: side!(|x, _| x.exp()),
                libm_side: side!(|x, _| libm::exp(x)),
            },
            vec![
                ArgumentRange {
                    label: "log-probabilities, uniform in [-700, 0]",
                    first: |x| -700.0 * x,
                    second: unused,
                },
                ArgumentRange {
                    label: "log-ratios, uniform in [-20, 20]",
                    first: |x| -20.0 + 40.0 * x,
                    second: unused,
                },
                ArgumentRange {
                    label: "log-probabilities past underflow, uniform in [-745, -700]",
                    first: |x| -745.0 + 45.0 * x,
                    second: unused,
                },
            ],
        ),
        (
            Function {
                name: "powf",
                std_side: side!(|x, y| x.powf(y)),
                libm_side: side!(|x, y| libm::pow(x, y)),
            },
            vec![
                ArgumentRange {
                    label: "10^(-q/10), q uniform in [0, 93]",
                    first: |_| 10.0,
                    second: |x| -9.3 * x,
                },
                ArgumentRange {
                    label: "base log-uniform in [1e-30, 1], exponent uniform in [0.1, 10]",
                    first: |x| log_uniform(x, 1e-30, 1.0),
                    second: |x| 0.1 + 9.9 * ((x * 7919.0).fract()),
                },
            ],
        ),
        (
            Function {
                name: "log10",
                std_side: side!(|x, _| x.log10()),
                libm_side: side!(|x, _| libm::log10(x)),
            },
            vec![ArgumentRange {
                label: "probabilities, log-uniform in [1e-300, 1]",
                first: |x| log_uniform(x, 1e-300, 1.0),
                second: unused,
            }],
        ),
        (
            Function {
                name: "ln_1p",
                std_side: side!(|x, _| x.ln_1p()),
                libm_side: side!(|x, _| libm::log1p(x)),
            },
            vec![ArgumentRange {
                label: "e^(lo - hi) in a two-term log-sum-exp, log-uniform in [1.4e-4, 1]",
                first: |x| log_uniform(x, 1.4e-4, 1.0),
                second: unused,
            }],
        ),
        (
            // No shipped call uses `exp_m1`; it is measured because the plan's module would
            // offer it.
            Function {
                name: "exp_m1",
                std_side: side!(|x, _| x.exp_m1()),
                libm_side: side!(|x, _| libm::expm1(x)),
            },
            vec![ArgumentRange {
                label: "-d, d log-uniform in [1e-12, 50]",
                first: |x| -log_uniform(x, 1e-12, 50.0),
                second: unused,
            }],
        ),
        (
            Function {
                name: "powi",
                std_side: side!(|x, y| x.powi(y as i32)),
                libm_side: side!(|x, y| powi_loop(x, y as i32)),
            },
            vec![ArgumentRange {
                label: "base uniform in (0, 1), exponent uniform in 0..=300 (std vs a Rust loop)",
                first: |x| 1e-6 + (1.0 - 2e-6) * x,
                second: |x| ((x * 7919.0).fract() * 301.0).floor(),
            }],
        ),
        (
            Function {
                name: "sin",
                std_side: side!(|x, _| x.sin()),
                libm_side: side!(|x, _| libm::sin(x)),
            },
            vec![ArgumentRange {
                label: "pi x for a concentration x in [1e-3, 0.5): uniform in [3.1e-3, pi/2)",
                first: |x| 3.1e-3 + (std::f64::consts::FRAC_PI_2 - 3.1e-3) * x,
                second: unused,
            }],
        ),
    ]
}

/// The arguments of one range, as two columns.
fn arguments(range: &ArgumentRange) -> (Vec<f64>, Vec<f64>) {
    (0..ARGUMENTS_PER_RANGE)
        .map(|index| {
            // Evenly spaced points in [0, 1), visited in a scrambled order so neither the branch
            // predictor nor the cache sees a sweep. Multiplying by an odd number modulo 2^19 and
            // then reversing the 19 bits is a permutation of the indices; the reversal is what
            // breaks up the regular strides the multiplication alone leaves in the low bits.
            let bits = ARGUMENTS_PER_RANGE.trailing_zeros();
            let scrambled = ((index as u32).wrapping_mul(0x9E37_79B9)
                & (ARGUMENTS_PER_RANGE as u32 - 1))
                .reverse_bits()
                >> (u32::BITS - bits);
            let spaced = f64::from(scrambled) / ARGUMENTS_PER_RANGE as f64;
            ((range.first)(spaced), (range.second)(spaced))
        })
        .unzip()
}

/// Time one pass of `side` over the arguments, writing its outputs. Returns nanoseconds per call.
fn timed_pass(side: Pass, first: &[f64], second: &[f64], out: &mut [f64]) -> f64 {
    let start = Instant::now();
    side(first, second, out);
    black_box(&out);
    start.elapsed().as_nanos() as f64 / first.len() as f64
}

/// Median, fastest and slowest of a set of timings.
fn summary(mut timings: Vec<f64>) -> (f64, f64, f64) {
    timings.sort_by(f64::total_cmp);
    (
        timings[timings.len() / 2],
        timings[0],
        timings[timings.len() - 1],
    )
}

/// Size of the difference between two outputs, in units in the last place of the smaller one: two
/// adjacent values are one unit apart, including either side of a power of two, where the gap
/// above the larger would be twice the gap between them. Identical values (including two NaNs of
/// the same bits) are zero apart; a finite value against a non-finite one is infinitely apart.
fn ulps_apart(a: f64, b: f64) -> f64 {
    if a.to_bits() == b.to_bits() {
        return 0.0;
    }
    if !(a.is_finite() && b.is_finite()) {
        return f64::INFINITY;
    }
    let smaller = a.abs().min(b.abs());
    let ulp = f64::from_bits(smaller.to_bits() + 1) - smaller;
    (a - b).abs() / ulp
}

/// The largest [`ulps_apart`] over two output columns, and how many pairs differ at all.
fn differences(a: &[f64], b: &[f64]) -> (usize, f64) {
    let differing = a
        .iter()
        .zip(b)
        .filter(|(x, y)| x.to_bits() != y.to_bits())
        .count();
    let largest = a
        .iter()
        .zip(b)
        .map(|(x, y)| ulps_apart(*x, *y))
        .fold(0.0, f64::max);
    (differing, largest)
}

fn file_name(function: &str, range_index: usize, side: &str) -> String {
    let function: String = function
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect();
    format!("{function}.{range_index}.{side}.f64")
}

fn write_outputs(path: &Path, values: &[f64]) {
    let bytes: Vec<u8> = values
        .iter()
        .flat_map(|value| value.to_le_bytes())
        .collect();
    std::fs::write(path, bytes).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
}

fn read_outputs(path: &Path) -> Vec<u64> {
    let bytes = std::fs::read(path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    let (whole, rest) = bytes.as_chunks::<8>();
    assert!(
        rest.is_empty(),
        "{}: not a whole number of f64s",
        path.display()
    );
    whole
        .iter()
        .map(|chunk| u64::from_le_bytes(*chunk))
        .collect()
}

fn measure(output_dir: &Path, trials: usize) {
    std::fs::create_dir_all(output_dir).expect("create the output directory");
    println!(
        "platform: {} {}; {ARGUMENTS_PER_RANGE} arguments a range; {trials} timed passes a side",
        std::env::consts::OS,
        std::env::consts::ARCH
    );
    println!(
        "function\trange\tstd ns/call median [min, max]\tlibm ns/call median [min, max]\t\
         libm/std\targuments differing\tlargest difference (ulps)"
    );
    for (function, ranges) in cases() {
        for (range_index, range) in ranges.iter().enumerate() {
            let (first, second) = arguments(range);
            let mut std_out = vec![0.0; first.len()];
            let mut libm_out = vec![0.0; first.len()];
            // Warm both sides, untimed, so page faults, lazy binding and the processor's climb
            // out of its idle state land outside the timings. One pass was not enough: the first
            // range measured after a quiet spell came out a third slower than in later runs.
            for _ in 0..WARM_UP_PASSES {
                timed_pass(function.std_side, &first, &second, &mut std_out);
                timed_pass(function.libm_side, &first, &second, &mut libm_out);
            }
            let mut std_timings = Vec::with_capacity(trials);
            let mut libm_timings = Vec::with_capacity(trials);
            for trial in 0..trials {
                // Alternate which side goes first, so a drift in machine speed hits both alike.
                if trial % 2 == 0 {
                    std_timings.push(timed_pass(function.std_side, &first, &second, &mut std_out));
                    libm_timings.push(timed_pass(
                        function.libm_side,
                        &first,
                        &second,
                        &mut libm_out,
                    ));
                } else {
                    libm_timings.push(timed_pass(
                        function.libm_side,
                        &first,
                        &second,
                        &mut libm_out,
                    ));
                    std_timings.push(timed_pass(function.std_side, &first, &second, &mut std_out));
                }
            }
            let (differing, largest) = differences(&std_out, &libm_out);
            let (std_median, std_min, std_max) = summary(std_timings);
            let (libm_median, libm_min, libm_max) = summary(libm_timings);
            println!(
                "{}\t{}\t{std_median:.2} [{std_min:.2}, {std_max:.2}]\t\
                 {libm_median:.2} [{libm_min:.2}, {libm_max:.2}]\t{:.2}\t{differing}\t{largest}",
                function.name,
                range.label,
                libm_median / std_median,
            );
            write_outputs(
                &output_dir.join(file_name(function.name, range_index, "std")),
                &std_out,
            );
            write_outputs(
                &output_dir.join(file_name(function.name, range_index, "libm")),
                &libm_out,
            );
        }
    }
    check_folded_calls(output_dir);
}

/// Calls written with literal operands, so the compiler may compute them while building. LLVM's
/// constant folder evaluates `ln`, `exp`, `pow` and `sin` of constants with the maths library of
/// the machine running the compiler, and `powi` with that machine's `pow`, so a folded value can
/// differ from the same call made at run time — and from one build machine to another. Each entry
/// is (name, folded at build time, std at run time, libm or the Rust loop at run time).
fn folded_calls() -> Vec<(String, f64, f64, f64)> {
    macro_rules! unary {
        ($method:ident, $libm:path; $($x:literal),* $(,)?) => {
            vec![$((
                format!("{}({})", stringify!($method), $x),
                f64::$method($x),
                black_box($x as f64).$method(),
                $libm(black_box($x)),
            )),*]
        };
    }
    macro_rules! powf {
        ($(($base:literal, $exponent:literal)),* $(,)?) => {
            vec![$((
                format!("powf({}, {})", $base, $exponent),
                f64::powf($base, $exponent),
                black_box($base as f64).powf(black_box($exponent)),
                libm::pow(black_box($base), black_box($exponent)),
            )),*]
        };
    }
    macro_rules! powi {
        ($(($base:literal, $exponent:literal)),* $(,)?) => {
            vec![$((
                format!("powi({}, {})", $base, $exponent),
                f64::powi($base, $exponent),
                black_box($base as f64).powi(black_box($exponent)),
                powi_loop(black_box($base), black_box($exponent)),
            )),*]
        };
    }
    let mut calls = Vec::new();
    // The literal arguments of the caller's constant sites (inventory §2: the shipped stutter
    // shares, ln 3, the quality ladder's base) and some others.
    calls.extend(unary!(ln, libm::log;
        0.88, 0.0475, 0.05, 3.0, 2.0, 10.0, 0.3, 0.123_456_789, 1e-12, 1e-300, 0.999, 7.5,
        1234.5678, 0.01, 0.99, 0.2, 0.8, 1.0e8, 0.6932, 42.0));
    calls.extend(unary!(exp, libm::exp;
        -0.5, -1.0, -2.45, -10.0, -20.5, -100.25, -700.0, 0.3, 1.5, 3.7, -0.001, -5.55,
        -37.2, -250.0, 2.0, -0.123_456_789, -1e-6, -60.0, -3.3, -13.0));
    calls.extend(unary!(sin, libm::sin;
        0.001, 0.01, 0.1, 0.5, 1.0, 1.5, 0.314_159, 0.77, 1.2, 0.0031, 0.25, 0.75, 1.05,
        0.42, 0.9, 1.3, 0.066, 0.2, 0.6, 1.57));
    calls.extend(powf!(
        (10.0, -0.4),
        (10.0, -1.3),
        (10.0, -2.7),
        (10.0, -9.3),
        (1.47, 0.5),
        (0.999, 150.5),
        (0.95, 3.3),
        (0.5, 0.25),
        (2.0, 0.1),
        (0.123_456_789, 2.5),
        (7.0, 0.333_333_333_333),
        (0.001, 0.7),
        (1.15, 23.5),
        (0.6, 7.2),
        (0.05, 9.1),
        (0.2, 11.9),
        (3.0, 0.31),
        (0.875, 6.3),
        (0.1, 3.0),
        (1.000_001, 250.0),
    ));
    calls.extend(powi!(
        (0.5, 3),
        (0.5, 17),
        (0.95, 37),
        (0.999, 150),
        (0.999, 299),
        (0.99, 5),
        (0.55, 4),
        (0.6, 7),
        (0.05, 9),
        (0.2, 11),
        (1.47, 12),
        (1.15, 23),
        (0.9, 101),
        (0.123_456_789, 13),
        (0.7, -9),
        (1.000_001, 250),
        (0.333_333_333_333, 29),
        (0.875, 63),
        (0.1, 30),
        (3.0, 31),
    ));
    calls
}

/// Report, per function, how many literal calls folded at build time differ from the same call at
/// run time through std, and from libm (or, for `powi`, the Rust loop); print every one that
/// differs from std's run-time value; and write the folded values so `compare` can set two
/// builds side by side.
fn check_folded_calls(output_dir: &Path) {
    let calls = folded_calls();
    println!(
        "function\tliteral calls\tfolded != std at run time\tfolded != libm or loop\tlargest (ulps)"
    );
    for function in ["ln", "exp", "sin", "powf", "powi"] {
        let of_function: Vec<_> = calls
            .iter()
            .filter(|(name, ..)| name.starts_with(&format!("{function}(")))
            .collect();
        let against_std = of_function
            .iter()
            .filter(|(_, folded, run, _)| folded.to_bits() != run.to_bits())
            .count();
        let against_libm = of_function
            .iter()
            .filter(|(_, folded, _, libm)| folded.to_bits() != libm.to_bits())
            .count();
        let largest = of_function
            .iter()
            .map(|(_, folded, run, _)| ulps_apart(*folded, *run))
            .fold(0.0, f64::max);
        println!(
            "{function}\t{}\t{against_std}\t{against_libm}\t{largest}",
            of_function.len()
        );
    }
    for (name, folded, run, libm) in &calls {
        if folded.to_bits() != run.to_bits() {
            println!(
                "  {name}: folded {:#018x}, std at run time {:#018x}, libm or loop {:#018x}",
                folded.to_bits(),
                run.to_bits(),
                libm.to_bits()
            );
        }
    }
    let values: Vec<f64> = calls.iter().map(|(_, folded, ..)| *folded).collect();
    write_outputs(&output_dir.join("folded_calls.f64"), &values);
}

fn compare(dir_a: &Path, dir_b: &Path) {
    let folded_a = read_outputs(&dir_a.join("folded_calls.f64"));
    let folded_b = read_outputs(&dir_b.join("folded_calls.f64"));
    println!(
        "literal calls folded at build time: {} of {} differ between the two builds",
        folded_a
            .iter()
            .zip(&folded_b)
            .filter(|(a, b)| a != b)
            .count(),
        folded_a.len()
    );
    println!(
        "function\trange\tstd outputs differing between the two\tlargest (ulps)\t\
         libm outputs differing\tlargest (ulps)"
    );
    for (function, ranges) in cases() {
        for (range_index, range) in ranges.iter().enumerate() {
            let side = |side: &str| {
                let name = file_name(function.name, range_index, side);
                let a: Vec<f64> = read_outputs(&dir_a.join(&name))
                    .into_iter()
                    .map(f64::from_bits)
                    .collect();
                let b: Vec<f64> = read_outputs(&dir_b.join(&name))
                    .into_iter()
                    .map(f64::from_bits)
                    .collect();
                assert_eq!(
                    a.len(),
                    b.len(),
                    "{name}: the two runs used different argument counts"
                );
                differences(&a, &b)
            };
            let (std_differing, std_largest) = side("std");
            let (libm_differing, libm_largest) = side("libm");
            println!(
                "{}\t{}\t{std_differing}\t{std_largest}\t{libm_differing}\t{libm_largest}",
                function.name, range.label
            );
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("measure") if args.len() >= 3 => {
            let trials = args
                .get(3)
                .map_or(DEFAULT_TRIALS, |t| t.parse().expect("trials"));
            measure(Path::new(&args[2]), trials);
        }
        Some("compare") if args.len() == 4 => compare(Path::new(&args[2]), Path::new(&args[3])),
        _ => {
            eprintln!(
                "usage: float_libm_vs_std measure <output-dir> [trials]\n\
                 \x20      float_libm_vs_std compare <dir-a> <dir-b>"
            );
            std::process::exit(2);
        }
    }
}
