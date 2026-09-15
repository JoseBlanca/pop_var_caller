// Ported from musl's src/math/exp.c and src/math/exp_data.c, which are
// Copyright (c) 2018, Arm Limited. SPDX-License-Identifier: MIT
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software
// and associated documentation files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or
// substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
// BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

//! `eˣ` by table lookup, the algorithm glibc and musl use, written in Rust.
//!
//! **Why this is not `libm::exp`.** The `libm` crate's `exp` is FreeBSD's rational-function
//! version. On arm64 Linux it takes 1.3 to 2.1 times as long a call as glibc's, except past
//! underflow, where the two are equal (report `portable_float_A2_libm_vs_std_2026-09-14.md`).
//! `exp` took about a fifth of the parameter fit's CPU time before the conversion, the most of any
//! maths function (report `portable_float_A3_caller_baseline_2026-09-14.md`). glibc (since 2.28)
//! and musl compute `exp` with Arm's table-driven algorithm (optimized-routines, 2018, MIT); this
//! module is that algorithm, ported line for line from musl's `src/math/exp.c`. Being plain Rust
//! arithmetic, it gives the same bits on every platform. glibc's build of the same algorithm does
//! not match it exactly: on arm64 Linux it gave different bits on 274 to 404 of 524,288 arguments
//! in each of three ranges (`tmp/exp_D1/linux`). Arm's source is written to be built with or without
//! fused multiply-add, which Rust never uses; that is the likely reason, not a measured one.
//!
//! **How it works.** Write `x = k·ln2/128 + r` with `k` an integer and `|r| ≤ ln2/256`. Then
//! `eˣ = 2^(k/128) · eʳ`. The first factor is `2^(k/128 mod 1)` from a 128-entry table, scaled by a
//! power of two, which is exact; the second is a degree-5 polynomial in `r`, accurate because `r`
//! is tiny. Each table entry is held as the nearest double to `2^(j/128)` plus that double's
//! rounding error, relative and itself rounded, so the lookup is exact to about 2⁻¹⁰⁶. musl's data
//! file bounds the error at 0.511 of a unit in the last place without fused multiply-add
//! (`exp_data.c`: "ulp error: 0.509 (0.511 without fma)"); `libm`'s own `exp` is bounded at under
//! one unit. Measured against correctly rounded values over 32,768 arguments in each of three wide
//! ranges, this port rounds 26 to 40 the other way, `libm` 3,152 to 3,220, never by more than one
//! step either (review `portable_float_D1_review_2026-09-15.md`).
//!
//! **The table is generated, not copied**: `scripts/float_exp_table.py` computes it at 200 bits,
//! checks six of its 128 entries against musl's `exp_data.c`, and with `--check` confirms the table
//! below is its output.

/// Bits of the table index: 128 entries.
const TABLE_BITS: u32 = 7;
/// Entries in the table.
const N: u64 = 1 << TABLE_BITS;
/// `128 / ln 2`, which is `0x1.71547652b82fep7`.
const INV_LN2_N: f64 = f64::from_bits(0x4067_1547_652b_82fe);
/// `−ln 2 / 128`, high part, `-0x1.62e42fefa0000p-8`: its low bits are zero so that `k` times it is
/// exact for every `k` this function meets.
const NEG_LN2_HI_N: f64 = f64::from_bits(0xbf76_2e42_fefa_0000);
/// `−ln 2 / 128`, low part, `-0x1.cf79abc9e3b3ap-47`.
const NEG_LN2_LO_N: f64 = f64::from_bits(0xbd0c_f79a_bc9e_3b3a);
/// `1.5 · 2⁵²`, `0x1.8p52`: added before rounding so that `z + SHIFT` is a positive double whose
/// low 52 bits hold `2⁵¹ + round(z)`; only `round(z) mod 128` and those bits shifted into `top` are
/// used.
const SHIFT: f64 = f64::from_bits(0x4338_0000_0000_0000);
/// The polynomial for `eʳ − 1 − r`, from `r²` up: `0x1.ffffffffffdbdp-2`, `0x1.555555555543cp-3`,
/// `0x1.55555cf172b91p-5`, `0x1.1111167a4d017p-7`.
const C2: f64 = f64::from_bits(0x3fdf_ffff_ffff_fdbd);
const C3: f64 = f64::from_bits(0x3fc5_5555_5555_543c);
const C4: f64 = f64::from_bits(0x3fa5_5555_cf17_2b91);
const C5: f64 = f64::from_bits(0x3f81_1111_67a4_d017);
/// `2¹⁰⁰⁹` and `2⁻¹⁰²²`, for results near overflow and underflow.
const TWO_TO_1009: f64 = f64::from_bits(0x7f00_0000_0000_0000);
const TWO_TO_MINUS_1022: f64 = f64::from_bits(0x0010_0000_0000_0000);

/// The top 12 bits of a double: its sign and exponent.
#[inline]
fn top12(x: f64) -> u32 {
    (x.to_bits() >> 52) as u32
}

/// `e` raised to `x`.
#[inline]
pub(super) fn exp(x: f64) -> f64 {
    // `0x3c9` is the exponent of 2⁻⁵⁴, `0x408` of 512 and `0x409` of 1024.
    let mut abstop = top12(x) & 0x7ff;
    if abstop.wrapping_sub(0x3c9) >= 0x408 - 0x3c9 {
        if abstop.wrapping_sub(0x3c9) >= 0x8000_0000 {
            // |x| < 2⁻⁵⁴, including zero: 1 + x rounds to the right answer.
            return 1.0 + x;
        }
        if abstop >= 0x409 {
            if x.to_bits() == f64::NEG_INFINITY.to_bits() {
                return 0.0;
            }
            if abstop >= 0x7ff {
                // +∞ or NaN.
                return 1.0 + x;
            }
            return if x.to_bits() >> 63 != 0 {
                0.0
            } else {
                f64::INFINITY
            };
        }
        // 512 ≤ |x| < 1024: the result may overflow or underflow, handled below.
        abstop = 0;
    }

    // x = k·ln2/N + r, with integer k and r in [−ln2/2N, ln2/2N].
    let z = INV_LN2_N * x;
    let mut kd = z + SHIFT;
    let ki = kd.to_bits();
    kd -= SHIFT;
    let r = x + kd * NEG_LN2_HI_N + kd * NEG_LN2_LO_N;
    // 2^(k/N) ≈ scale · (1 + tail).
    let index = (2 * (ki % N)) as usize;
    let top = ki << (52 - TABLE_BITS);
    let tail = f64::from_bits(TABLE[index]);
    // A valid scale only while −1023·N < k < 1024·N, which holds for |x| < 512; larger |x| goes
    // to `near_the_limits`, which moves the exponent back into range before using it.
    let scale_bits = TABLE[index + 1].wrapping_add(top);
    // eˣ = 2^(k/N) · eʳ ≈ scale + scale · (tail + eʳ − 1).
    let r2 = r * r;
    let tmp = tail + r + r2 * (C2 + r * C3) + r2 * r2 * (C4 + r * C5);
    if abstop == 0 {
        return near_the_limits(tmp, scale_bits, ki);
    }
    let scale = f64::from_bits(scale_bits);
    scale + scale * tmp
}

/// `scale · (1 + tmp)` when the scale's exponent may have left the double's range: `k > 0` may
/// overflow, `k < 0` may underflow into the subnormals, where the sum is rounded once, carefully,
/// before scaling, so the result is not rounded twice.
#[cold]
fn near_the_limits(tmp: f64, scale_bits: u64, ki: u64) -> f64 {
    if ki & 0x8000_0000 == 0 {
        // k > 0: the scale's exponent may have overflowed by up to 460.
        let scale = f64::from_bits(scale_bits.wrapping_sub(1009 << 52));
        return TWO_TO_1009 * (scale + scale * tmp);
    }
    // k < 0.
    let scale = f64::from_bits(scale_bits.wrapping_add(1022 << 52));
    let mut y = scale + scale * tmp;
    if y < 1.0 {
        let mut lo = scale - y + scale * tmp;
        let hi = 1.0 + y;
        // musl writes `lo = 1.0 - hi + y + lo`; adding `lo` last from either side is the same
        // IEEE addition, which is commutative.
        lo += 1.0 - hi + y;
        y = (hi + lo) - 1.0;
        // No −0.0: the sum above can round to it only under directed rounding, but keep the
        // result positive whatever produced it.
        if y == 0.0 {
            y = 0.0;
        }
    }
    TWO_TO_MINUS_1022 * y
}

/// `2^(j/128)` for `j` in `0..128`, as pairs: the bits of the rounding error of the nearest double,
/// relative to it, and the bits of that nearest double minus `j·2⁴⁵`, so that adding `k·2⁴⁵`
/// (wrapping) gives the bits of `2^(k/128)`'s scale for any integer `−1023·128 < k < 1024·128`.
// Generated by scripts/float_exp_table.py; do not edit by hand.
#[rustfmt::skip]
const TABLE: [u64; 2 * N as usize] = [
    0x0000000000000000, 0x3ff0000000000000,
    0x3c9b3b4f1a88bf6e, 0x3feff63da9fb3335,
    0xbc7160139cd8dc5d, 0x3fefec9a3e778061,
    0xbc905e7a108766d1, 0x3fefe315e86e7f85,
    0x3c8cd2523567f613, 0x3fefd9b0d3158574,
    0xbc8bce8023f98efa, 0x3fefd06b29ddf6de,
    0x3c60f74e61e6c861, 0x3fefc74518759bc8,
    0x3c90a3e45b33d399, 0x3fefbe3ecac6f383,
    0x3c979aa65d837b6d, 0x3fefb5586cf9890f,
    0x3c8eb51a92fdeffc, 0x3fefac922b7247f7,
    0x3c3ebe3d702f9cd1, 0x3fefa3ec32d3d1a2,
    0xbc6a033489906e0b, 0x3fef9b66affed31b,
    0xbc9556522a2fbd0e, 0x3fef9301d0125b51,
    0xbc5080ef8c4eea55, 0x3fef8abdc06c31cc,
    0xbc91c923b9d5f416, 0x3fef829aaea92de0,
    0x3c80d3e3e95c55af, 0x3fef7a98c8a58e51,
    0xbc801b15eaa59348, 0x3fef72b83c7d517b,
    0xbc8f1ff055de323d, 0x3fef6af9388c8dea,
    0x3c8b898c3f1353bf, 0x3fef635beb6fcb75,
    0xbc96d99c7611eb26, 0x3fef5be084045cd4,
    0x3c9aecf73e3a2f60, 0x3fef54873168b9aa,
    0xbc8fe782cb86389d, 0x3fef4d5022fcd91d,
    0x3c8a6f4144a6c38d, 0x3fef463b88628cd6,
    0x3c807a05b0e4047d, 0x3fef3f49917ddc96,
    0x3c968efde3a8a894, 0x3fef387a6e756238,
    0x3c875e18f274487d, 0x3fef31ce4fb2a63f,
    0x3c80472b981fe7f2, 0x3fef2b4565e27cdd,
    0xbc96b87b3f71085e, 0x3fef24dfe1f56381,
    0x3c82f7e16d09ab31, 0x3fef1e9df51fdee1,
    0xbc3d219b1a6fbffa, 0x3fef187fd0dad990,
    0x3c8b3782720c0ab4, 0x3fef1285a6e4030b,
    0x3c6e149289cecb8f, 0x3fef0cafa93e2f56,
    0x3c834d754db0abb6, 0x3fef06fe0a31b715,
    0x3c864201e2ac744c, 0x3fef0170fc4cd831,
    0x3c8fdd395dd3f84a, 0x3feefc08b26416ff,
    0xbc86a3803b8e5b04, 0x3feef6c55f929ff1,
    0xbc924aedcc4b5068, 0x3feef1a7373aa9cb,
    0xbc9907f81b512d8e, 0x3feeecae6d05d866,
    0xbc71d1e83e9436d2, 0x3feee7db34e59ff7,
    0xbc991919b3ce1b15, 0x3feee32dc313a8e5,
    0x3c859f48a72a4c6d, 0x3feedea64c123422,
    0xbc9312607a28698a, 0x3feeda4504ac801c,
    0xbc58a78f4817895b, 0x3feed60a21f72e2a,
    0xbc7c2c9b67499a1b, 0x3feed1f5d950a897,
    0x3c4363ed60c2ac11, 0x3feece086061892d,
    0x3c9666093b0664ef, 0x3feeca41ed1d0057,
    0x3c6ecce1daa10379, 0x3feec6a2b5c13cd0,
    0x3c93ff8e3f0f1230, 0x3feec32af0d7d3de,
    0x3c7690cebb7aafb0, 0x3feebfdad5362a27,
    0x3c931dbdeb54e077, 0x3feebcb299fddd0d,
    0xbc8f94340071a38e, 0x3feeb9b2769d2ca7,
    0xbc87deccdc93a349, 0x3feeb6daa2cf6642,
    0xbc78dec6bd0f385f, 0x3feeb42b569d4f82,
    0xbc861246ec7b5cf6, 0x3feeb1a4ca5d920f,
    0x3c93350518fdd78e, 0x3feeaf4736b527da,
    0x3c7b98b72f8a9b05, 0x3feead12d497c7fd,
    0x3c9063e1e21c5409, 0x3feeab07dd485429,
    0x3c34c7855019c6ea, 0x3feea9268a5946b7,
    0x3c9432e62b64c035, 0x3feea76f15ad2148,
    0xbc8ce44a6199769f, 0x3feea5e1b976dc09,
    0xbc8c33c53bef4da8, 0x3feea47eb03a5585,
    0xbc845378892be9ae, 0x3feea34634ccc320,
    0xbc93cedd78565858, 0x3feea23882552225,
    0x3c5710aa807e1964, 0x3feea155d44ca973,
    0xbc93b3efbf5e2228, 0x3feea09e667f3bcd,
    0xbc6a12ad8734b982, 0x3feea012750bdabf,
    0xbc6367efb86da9ee, 0x3fee9fb23c651a2f,
    0xbc80dc3d54e08851, 0x3fee9f7df9519484,
    0xbc781f647e5a3ecf, 0x3fee9f75e8ec5f74,
    0xbc86ee4ac08b7db0, 0x3fee9f9a48a58174,
    0xbc8619321e55e68a, 0x3fee9feb564267c9,
    0x3c909ccb5e09d4d3, 0x3feea0694fde5d3f,
    0xbc7b32dcb94da51d, 0x3feea11473eb0187,
    0x3c94ecfd5467c06b, 0x3feea1ed0130c132,
    0x3c65ebe1abd66c55, 0x3feea2f336cf4e62,
    0xbc88a1c52fb3cf42, 0x3feea427543e1a12,
    0xbc9369b6f13b3734, 0x3feea589994cce13,
    0xbc805e843a19ff1e, 0x3feea71a4623c7ad,
    0xbc94d450d872576e, 0x3feea8d99b4492ed,
    0x3c90ad675b0e8a00, 0x3feeaac7d98a6699,
    0x3c8db72fc1f0eab4, 0x3feeace5422aa0db,
    0xbc65b6609cc5e7ff, 0x3feeaf3216b5448c,
    0x3c7bf68359f35f44, 0x3feeb1ae99157736,
    0xbc93091fa71e3d83, 0x3feeb45b0b91ffc6,
    0xbc5da9b88b6c1e29, 0x3feeb737b0cdc5e5,
    0xbc6c23f97c90b959, 0x3feeba44cbc8520f,
    0xbc92434322f4f9aa, 0x3feebd829fde4e50,
    0xbc85ca6cd7668e4b, 0x3feec0f170ca07ba,
    0x3c71affc2b91ce27, 0x3feec49182a3f090,
    0x3c6dd235e10a73bb, 0x3feec86319e32323,
    0xbc87c50422622263, 0x3feecc667b5de565,
    0x3c8b1c86e3e231d5, 0x3feed09bec4a2d33,
    0xbc91bbd1d3bcbb15, 0x3feed503b23e255d,
    0x3c90cc319cee31d2, 0x3feed99e1330b358,
    0x3c8469846e735ab3, 0x3feede6b5579fdbf,
    0xbc82dfcd978e9db4, 0x3feee36bbfd3f37a,
    0x3c8c1a7792cb3387, 0x3feee89f995ad3ad,
    0xbc907b8f4ad1d9fa, 0x3feeee07298db666,
    0xbc55c3d956dcaeba, 0x3feef3a2b84f15fb,
    0xbc90a40e3da6f640, 0x3feef9728de5593a,
    0xbc68d6f438ad9334, 0x3feeff76f2fb5e47,
    0xbc91eee26b588a35, 0x3fef05b030a1064a,
    0x3c74ffd70a5fddcd, 0x3fef0c1e904bc1d2,
    0xbc91bdfbfa9298ac, 0x3fef12c25bd71e09,
    0x3c736eae30af0cb3, 0x3fef199bdd85529c,
    0x3c8ee3325c9ffd94, 0x3fef20ab5fffd07a,
    0x3c84e08fd10959ac, 0x3fef27f12e57d14b,
    0x3c63cdaf384e1a67, 0x3fef2f6d9406e7b5,
    0x3c676b2c6c921968, 0x3fef3720dcef9069,
    0xbc808a1883ccb5d2, 0x3fef3f0b555dc3fa,
    0xbc8fad5d3ffffa6f, 0x3fef472d4a07897c,
    0xbc900dae3875a949, 0x3fef4f87080d89f2,
    0x3c74a385a63d07a7, 0x3fef5818dcfba487,
    0xbc82919e2040220f, 0x3fef60e316c98398,
    0x3c8e5a50d5c192ac, 0x3fef69e603db3285,
    0x3c843a59ac016b4b, 0x3fef7321f301b460,
    0xbc82d52107b43e1f, 0x3fef7c97337b9b5f,
    0xbc892ab93b470dc9, 0x3fef864614f5a129,
    0x3c74b604603a88d3, 0x3fef902ee78b3ff6,
    0x3c83c5ec519d7271, 0x3fef9a51fbc74c83,
    0xbc8ff7128fd391f0, 0x3fefa4afa2a490da,
    0xbc8dae98e223747d, 0x3fefaf482d8e67f1,
    0x3c8ec3bc41aa2008, 0x3fefba1bee615a27,
    0x3c842b94c3a9eb32, 0x3fefc52b376bba97,
    0x3c8a64a931d185ee, 0x3fefd0765b6e4540,
    0xbc8e37bae43be3ed, 0x3fefdbfdad9cbe14,
    0x3c77893b4d91cd9d, 0x3fefe7c1819e90d8,
    0x3c5305c14160cc89, 0x3feff3c22b8f71f1,
];

#[cfg(test)]
mod tests {
    use super::*;
    use std::hint::black_box;

    /// Steps between two doubles of one sign: adjacent doubles are one apart.
    fn steps_apart(a: f64, b: f64) -> u64 {
        a.to_bits().abs_diff(b.to_bits())
    }

    /// **Every table entry is what its doc says**: adding `j·2⁴⁵` to the high half gives a double
    /// within one step of `2^(j/128)`, and the tail is a relative correction smaller than half a
    /// step. Checked against `libm::pow`, which is within a step of the true value; a mistyped or
    /// misordered entry is off by far more.
    #[test]
    fn every_table_entry_is_two_to_the_j_over_128() {
        for j in 0..N {
            let scale = f64::from_bits(TABLE[(2 * j + 1) as usize].wrapping_add(j << 45));
            let exact = libm::pow(2.0, j as f64 / N as f64);
            assert!(
                steps_apart(scale, exact) <= 1,
                "entry {j}: {scale:e} against 2^({j}/128) = {exact:e}"
            );
            let tail = f64::from_bits(TABLE[(2 * j) as usize]);
            assert!(tail.abs() <= f64::EPSILON / 2.0, "entry {j}: tail {tail:e}");
        }
    }

    /// Within one step of `libm::exp` over the whole range where `eˣ` is a finite non-zero double,
    /// on a dense grid. Both are bounded under one unit in the last place of the true value, so a
    /// mistyped table entry or a wrong threshold fails here with a readable message. **An error
    /// smaller than a step does not** — a slightly wrong coefficient or tail —
    /// `the_grid_has_musls_bits` catches those.
    #[test]
    fn agrees_with_libm_to_within_one_step_everywhere_it_is_finite() {
        let (low, high) = (-745.13, 709.78);
        let points = 200_000;
        for i in 0..=points {
            let x = black_box(low + (high - low) * f64::from(i) / f64::from(points));
            let ours = exp(x);
            let theirs = libm::exp(x);
            assert!(
                steps_apart(ours, theirs) <= 1,
                "exp({x}): {ours:e} against libm's {theirs:e}"
            );
        }
    }

    /// **The same grid's bits, folded into one number.** The expected value was computed from an
    /// independent transliteration of musl's `exp.c` run on IEEE doubles (review
    /// `portable_float_D1_review_2026-09-15.md`, `tmp/review_D/verify_port.py`), so a coefficient,
    /// table tail or subnormal-rounding error that stays within one step still fails here. The
    /// review's mutations give, instead: `C4` set to 1/24 `0x9a73d4d8007071ed`, entry 37's tail
    /// zeroed `0xe59df08c2913e003`, the subnormal rounding removed `0x7757c8921cbcd118`. **If this
    /// fails, the port and musl differ somewhere; find it rather than re-recording.**
    #[test]
    fn the_grid_has_musls_bits() {
        let (low, high) = (-745.13, 709.78);
        let points = 200_000;
        let mut digest: u64 = 0xcbf2_9ce4_8422_2325;
        for i in 0..=points {
            let x = black_box(low + (high - low) * f64::from(i) / f64::from(points));
            digest = (digest ^ exp(x).to_bits()).wrapping_mul(0x0100_0000_01b3);
        }
        assert_eq!(digest, 0xb1d1_6a18_c6a4_2d12);
    }

    /// Arguments that reach each path's boundary, with bits checked against correct rounding at
    /// 300 bits by the review: the subnormal rounding step, one table entry's tail, subnormal
    /// arguments, the first arguments of the slow path, the last finite and first infinite result,
    /// the last non-zero and first zero result, and the |x| ≥ 1024 branch with a finite argument.
    #[test]
    fn keeps_the_bits_at_each_paths_boundary() {
        let cases: [(f64, u64); 10] = [
            (-711.252_420_65, 0x0000_eb83_3419_deff),
            (0.200_362_856_880_609_17, 0x3ff3_8cae_6d05_d865),
            (f64::from_bits(1), 1.0_f64.to_bits()),
            (-f64::from_bits(1), 1.0_f64.to_bits()),
            (512.0, 0x6e19_4765_04ba_852e),
            (-512.0, 0x11c4_4109_edb2_0931),
            (f64::from_bits(0x4086_2e42_fefa_39ef), 0x7fef_ffff_ffff_ff2a),
            (
                f64::from_bits(0x4086_2e42_fefa_39f0),
                f64::INFINITY.to_bits(),
            ),
            (f64::from_bits(0xc087_4910_d52d_3051), 0x0000_0000_0000_0001),
            (f64::from_bits(0xc087_4910_d52d_3052), 0.0_f64.to_bits()),
        ];
        for (x, bits) in cases {
            let ours = exp(black_box(x));
            assert_eq!(
                ours.to_bits(),
                bits,
                "exp({x:e}) = {ours:e}, expected bits {bits:#018x}"
            );
        }
        assert_eq!(exp(black_box(f64::MAX)), f64::INFINITY);
        assert_eq!(exp(black_box(-f64::MAX)).to_bits(), 0.0_f64.to_bits());
        assert!(exp(black_box(-f64::NAN)).is_nan());
    }

    /// The edges: zero of either sign, the arguments too small to move `1`, NaN, the infinities,
    /// the last finite results before overflow, and the subnormal tail before underflow.
    #[test]
    fn keeps_its_edges() {
        assert_eq!(exp(black_box(0.0)), 1.0);
        assert_eq!(exp(black_box(-0.0)), 1.0);
        assert_eq!(exp(black_box(1e-300)), 1.0);
        assert!(exp(black_box(f64::NAN)).is_nan());
        assert_eq!(exp(black_box(f64::INFINITY)), f64::INFINITY);
        assert_eq!(
            exp(black_box(f64::NEG_INFINITY)).to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(
            exp(black_box(709.78)).to_bits(),
            libm::exp(709.78).to_bits()
        );
        assert_eq!(exp(black_box(709.79)), f64::INFINITY);
        assert_eq!(exp(black_box(1000.0)), f64::INFINITY);
        assert_eq!(exp(black_box(-745.13)), f64::from_bits(1));
        assert_eq!(exp(black_box(-746.0)).to_bits(), 0.0_f64.to_bits());
        assert_eq!(exp(black_box(-1000.0)).to_bits(), 0.0_f64.to_bits());
        assert_eq!(exp(black_box(1.0)), std::f64::consts::E);
    }
}
