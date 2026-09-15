//! Floating-point functions that return the same bits on every platform.
//!
//! Rust's `f64::ln`, `exp`, `powf` and the other transcendental methods hand their work to the
//! operating system's maths library — glibc on Linux, Apple's libm on macOS — and IEEE 754 does not
//! require those to round alike. For `exp` they differ in the last binary place on between one
//! argument in 540 and one in 2,700, and a one-unit difference is enough to flip a tie: on 2026-09-13
//! one did, and a repeat tract measured a base longer on macOS than on Linux. A caller whose output
//! depends on the machine cannot be checked on one machine for another.
//!
//! Every function here computes in Rust instead — through the [`libm`] crate, for `exp` by a
//! table-driven algorithm written here, or for `powi` as plain multiplication — so the same argument
//! gives the same bits wherever the binary runs. On aarch64
//! and x86_64, libm's only processor-specific code is `sqrt`, `fma` and (on aarch64) `rint`, which
//! IEEE 754 requires to be rounded exactly; `pow` calls `sqrt`, and none of the functions this module
//! uses calls `fma`. The pinned bits in the tests were recorded on aarch64 macOS and aarch64 Linux; x86_64
//! is covered by the source, not yet by a run. The measurements behind this module — speed, and how
//! often each function's bits differ from each platform's library — are in
//! `doc/devel/reports/implementations/portable_float_A2_libm_vs_std_2026-09-14.md`.
//!
//! **What stays on std.** `+ − × ÷`, `sqrt`, `abs`, `floor` and the comparisons are rounded exactly
//! by IEEE 754 or involve no rounding, so std's are already portable. `ln Γ` is
//! [`genetics::lgamma`](crate::genetics::lgamma), which has always called libm.
//!
//! **What it costs.** libm's `ln` takes 0.5 to 1.5 ns a call more than the platform's, and `exp` up
//! to 3.3 ns more, on an Apple M5 Pro and in an arm64 Linux VM. That made the parameter fit about 30%
//! slower and the calling commands at most 3% slower (report
//! `portable_float_A3_caller_baseline_2026-09-14.md`), which is why `exp` was later replaced by the
//! table-driven version (`doc/devel/implementation_plans/portable_float.md`, Milestone D).
//!
//! **Nothing else may call std's versions.** `clippy.toml` refuses `f64::ln`, `exp`, `powf`, `powi`
//! and the other transcendental methods, and their `f32` twins, in every target, each with a message
//! pointing here. Two places allow them on purpose: this module's tests, which compare against
//! std, and the examples, which are research tools outside the guarantee.
//! `scripts/check_float_ban.sh` proves every entry still refuses its method.

mod table_exp;

/// The natural logarithm, `ln x`.
#[inline]
pub fn ln(x: f64) -> f64 {
    libm::log(x)
}

/// `e` raised to `x`.
///
/// **Computed by table lookup, not by `libm::exp`**: the algorithm glibc and musl use, written in
/// Rust in [`table_exp`]. `libm`'s rational-function `exp` takes 1.3 to 2.1 times as long a call as
/// glibc's outside underflow, and `exp` took about a fifth of the parameter fit's CPU time. Plain
/// Rust arithmetic, so the same argument gives the same bits on every platform.
#[inline]
pub fn exp(x: f64) -> f64 {
    table_exp::exp(x)
}

/// `base` raised to the real power `exponent`.
#[inline]
pub fn powf(base: f64, exponent: f64) -> f64 {
    libm::pow(base, exponent)
}

/// `base` raised to the integer power `exponent`, by repeated squaring.
///
/// **Written out rather than calling `f64::powi`.** At run time std's `powi` executes this same
/// square-and-multiply loop, and the two agreed on every argument measured, on both platforms. But
/// when both operands are known while building, the compiler computes std's `powi` itself, with the
/// build machine's `pow` — and that came out different from the run-time value in 15 of 20 cases, by
/// up to 26 units in the last place. Which calls get computed early depends on inlining, so the same
/// source line could give different bits in two builds. This loop is ordinary multiplication, which
/// the compiler can only evaluate exactly.
#[inline]
pub fn powi(base: f64, exponent: i32) -> f64 {
    let mut square = base;
    let mut remaining = exponent.unsigned_abs();
    let mut product = 1.0;
    loop {
        if remaining & 1 == 1 {
            product *= square;
        }
        remaining >>= 1;
        if remaining == 0 {
            break;
        }
        square *= square;
    }
    if exponent < 0 { 1.0 / product } else { product }
}

/// The base-10 logarithm, `log₁₀ x`.
#[inline]
pub fn log10(x: f64) -> f64 {
    libm::log10(x)
}

/// `ln(1 + x)`, accurate when `x` is near zero.
#[inline]
pub fn ln_1p(x: f64) -> f64 {
    libm::log1p(x)
}

/// `eˣ − 1`, accurate when `x` is near zero.
#[inline]
pub fn exp_m1(x: f64) -> f64 {
    libm::expm1(x)
}

/// The sine of `x` radians.
#[inline]
pub fn sin(x: f64) -> f64 {
    libm::sin(x)
}

/// The cosine of `x` radians.
#[inline]
pub fn cos(x: f64) -> f64 {
    libm::cos(x)
}

#[cfg(test)]
#[allow(
    clippy::disallowed_methods,
    reason = "these tests compare crate::float with std's methods on purpose"
)]
mod tests {
    use super::*;
    use std::hint::black_box;

    /// Outputs pinned as bits. They were produced by this module and were the same on both
    /// platforms recorded, so a change here means libm, the table-driven `exp`, or the function a
    /// name delegates to, changed.
    #[test]
    fn outputs_are_pinned_to_the_bit() {
        let pinned: [(&str, f64, u64); 16] = [
            ("ln(0.3)", ln(black_box(0.3)), PINNED[0]),
            ("ln(1e-300)", ln(black_box(1e-300)), PINNED[1]),
            ("ln(12345.678)", ln(black_box(12_345.678)), PINNED[2]),
            ("exp(-0.7)", exp(black_box(-0.7)), PINNED[3]),
            ("exp(-700.25)", exp(black_box(-700.25)), PINNED[4]),
            ("exp(13.5)", exp(black_box(13.5)), PINNED[5]),
            (
                "powf(10, -2.7)",
                powf(black_box(10.0), black_box(-2.7)),
                PINNED[6],
            ),
            (
                "powf(0.999, 150.5)",
                powf(black_box(0.999), black_box(150.5)),
                PINNED[7],
            ),
            (
                "powi(0.999, 299)",
                powi(black_box(0.999), black_box(299)),
                PINNED[8],
            ),
            (
                "powi(0.7, -9)",
                powi(black_box(0.7), black_box(-9)),
                PINNED[9],
            ),
            ("log10(3e-7)", log10(black_box(3e-7)), PINNED[10]),
            ("ln_1p(1.4e-4)", ln_1p(black_box(1.4e-4)), PINNED[11]),
            ("ln_1p(-0.3)", ln_1p(black_box(-0.3)), PINNED[12]),
            ("exp_m1(-1e-9)", exp_m1(black_box(-1e-9)), PINNED[13]),
            ("sin(0.77)", sin(black_box(0.77)), PINNED[14]),
            ("cos(1.2)", cos(black_box(1.2)), PINNED[15]),
        ];
        let moved: Vec<String> = pinned
            .iter()
            .filter(|(_, value, bits)| value.to_bits() != *bits)
            .map(|(name, value, bits)| {
                format!(
                    "{name} = {value:e}: bits {:#018x}, pinned {bits:#018x}",
                    value.to_bits()
                )
            })
            .collect();
        assert!(moved.is_empty(), "outputs moved:\n{}", moved.join("\n"));
    }

    /// Recorded on 2026-09-14 on macOS (aarch64) and in the Linux container (aarch64, glibc); both
    /// gave these bits. `exp(13.5)` was re-recorded on 2026-09-15 for the table-driven `exp`; the
    /// other two `exp` values did not move.
    const PINNED: [u64; 16] = [
        0xbff3_4378_fcbd_a721,
        0xc085_9634_47f8_7fb5,
        0x4022_d795_5979_1e31,
        0x3fdf_c80d_b9dd_5542,
        0x00ca_f5fe_9a48_5c8e,
        // exp(13.5): the table-driven `exp` gives the correctly rounded bits (0.44 of a step from
        // the true value, computed at 200 bits); libm's `exp` gave the next double up, `…ad8c`.
        0x4126_4290_bd5c_ad8b,
        0x3f60_585e_4c78_b079,
        0x3feb_86dd_50c8_818a,
        0x3fe7_b9f2_2a00_1e23,
        0x4038_c7eb_2c93_fcc8,
        0xc01a_176d_869b_02a0,
        0x3f22_594a_ab5c_4c32,
        0xbfd6_d3c3_24e1_3f4e,
        0xbe11_2e0b_e801_f1d9,
        0x3fe6_46bd_686f_ecc6,
        0x3fd7_30de_943b_79d4,
    ];

    /// Each function computes what its name says: within a few units in the last place of std's,
    /// which is the platform's library. Catches a name wired to the wrong libm function, which the
    /// pinned bits alone would only report as "changed".
    #[test]
    fn each_function_agrees_with_std_to_within_a_few_units_in_the_last_place() {
        /// Equal, or of one sign and at most four representable numbers apart. Counting steps
        /// rather than a relative tolerance keeps the slack at subnormal results, where
        /// `ε · |value|` rounds to zero. `+0.0 == −0.0` counts as equal.
        fn close(ours: f64, platform: f64) -> bool {
            if ours.to_bits() == platform.to_bits() || ours == platform {
                return true;
            }
            ours.is_sign_negative() == platform.is_sign_negative()
                && (ours.to_bits() as i64 - platform.to_bits() as i64).abs() <= 4
        }
        let x_values = black_box([
            1e-300, 1e-12, 0.001, 0.3, 0.999_999, 1.0, 2.0, 12_345.678, 1e12,
        ]);
        for x in x_values {
            assert!(close(ln(x), x.ln()), "ln({x})");
            assert!(close(log10(x), x.log10()), "log10({x})");
        }
        for x in black_box([-745.0, -700.25, -20.0, -1e-9, 0.0, 0.7, 13.5, 700.0]) {
            assert!(close(exp(x), x.exp()), "exp({x})");
            assert!(close(exp_m1(x), x.exp_m1()), "exp_m1({x})");
        }
        for x in black_box([-0.999, -0.3, -1e-12, 0.0, 1.4e-4, 1.0, 1e6]) {
            assert!(close(ln_1p(x), x.ln_1p()), "ln_1p({x})");
        }
        for x in black_box([0.0, 1e-3, 0.77, 1.2, std::f64::consts::FRAC_PI_2, 3.0, -2.5]) {
            assert!(close(sin(x), x.sin()), "sin({x})");
            assert!(close(cos(x), x.cos()), "cos({x})");
        }
        for (base, exponent) in black_box([
            (10.0, -2.7),
            (0.999, 150.5),
            (2.0, 0.5),
            (1e-30, 0.1),
            (7.0, 3.0),
        ]) {
            assert!(
                close(powf(base, exponent), base.powf(exponent)),
                "powf({base}, {exponent})"
            );
        }
    }

    /// The written-out `powi` gives std's run-time bits — the property that made writing it out
    /// change nothing at a call whose operands are not known while building.
    #[test]
    fn powi_matches_std_at_run_time_on_every_base_and_exponent_tried() {
        for base_index in 0..200 {
            let base = black_box(1e-6 + f64::from(base_index) * 0.012_345);
            for exponent in -40..=300 {
                let exponent = black_box(exponent);
                assert_eq!(
                    powi(base, exponent).to_bits(),
                    base.powi(exponent).to_bits(),
                    "powi({base}, {exponent})"
                );
            }
        }
    }

    /// The edges, against std at run time as well as against constants: extreme exponents, negative
    /// and signed-zero bases, infinities, NaN, and results that overflow before the reciprocal.
    #[test]
    fn powi_matches_std_at_run_time_on_its_edges() {
        let bases = [
            0.0,
            -0.0,
            -1.0,
            -0.5,
            -2.0,
            1.0,
            2.0,
            1e-200,
            1e200,
            f64::MIN_POSITIVE,
            f64::MAX,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN,
        ];
        let exponents = [
            i32::MIN,
            i32::MIN + 1,
            -1025,
            -2,
            -1,
            0,
            1,
            2,
            3,
            1024,
            i32::MAX,
        ];
        for base in bases {
            for exponent in exponents {
                let ours = powi(black_box(base), black_box(exponent));
                let platform = black_box(base).powi(black_box(exponent));
                assert!(
                    ours.to_bits() == platform.to_bits() || (ours.is_nan() && platform.is_nan()),
                    "powi({base}, {exponent}): {ours:e} against std's {platform:e}"
                );
            }
        }
    }

    #[test]
    fn powi_handles_its_edges() {
        assert_eq!(powi(0.37, 0), 1.0);
        assert_eq!(powi(0.0, 0), 1.0);
        assert_eq!(powi(2.0, 10), 1024.0);
        assert_eq!(powi(2.0, -2), 0.25);
        assert_eq!(powi(-2.0, 3), -8.0);
        assert_eq!(powi(0.0, -1), f64::INFINITY);
        assert_eq!(powi(2.0, i32::MIN), 0.0);
        assert!(powi(f64::NAN, 3).is_nan());
    }

    #[test]
    fn the_logarithms_and_exponentials_keep_their_limits() {
        assert_eq!(ln(0.0), f64::NEG_INFINITY);
        assert!(ln(-1.0).is_nan());
        assert_eq!(ln(1.0), 0.0);
        assert_eq!(exp(0.0), 1.0);
        assert_eq!(exp(f64::NEG_INFINITY), 0.0);
        assert_eq!(exp(710.0), f64::INFINITY);
        assert_eq!(ln_1p(-1.0), f64::NEG_INFINITY);
        assert_eq!(exp_m1(0.0), 0.0);
        assert_eq!(log10(1000.0), 3.0);
        assert_eq!(powf(2.0, 10.0), 1024.0);
    }
}
