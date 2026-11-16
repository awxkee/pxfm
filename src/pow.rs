/*
 * // Copyright (c) Radzivon Bartoshyk 4/2025. All rights reserved.
 * //
 * // Redistribution and use in source and binary forms, with or without modification,
 * // are permitted provided that the following conditions are met:
 * //
 * // 1.  Redistributions of source code must retain the above copyright notice, this
 * // list of conditions and the following disclaimer.
 * //
 * // 2.  Redistributions in binary form must reproduce the above copyright notice,
 * // this list of conditions and the following disclaimer in the documentation
 * // and/or other materials provided with the distribution.
 * //
 * // 3.  Neither the name of the copyright holder nor the names of its
 * // contributors may be used to endorse or promote products derived from
 * // this software without specific prior written permission.
 * //
 * // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * // AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * // IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * // DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * // FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * // DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * // SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * // CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * // OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
use crate::bits::{EXP_MASK, get_exponent_f64};
use crate::common::{f_fmla, is_integer};
use crate::double_double::DoubleDouble;
use crate::dyadic_float::{DyadicFloat128, DyadicSign};
use crate::exponents::exp;
use crate::logs::log_dyadic;
use crate::pow_exec::{exp_dyadic, pow_exp_1, pow_log_1};
use crate::triple_double::TripleDouble;
use crate::{f_exp2, f_exp10, log};

// The integer fast path only accepts negative exponents with hardware FMA.
const POWI_HAS_RECIP: bool = cfg!(any(
    all(
        any(target_arch = "x86", target_arch = "x86_64"),
        target_feature = "fma"
    ),
    target_arch = "aarch64"
));

#[inline]
fn is_odd_exponent(y_abs: u64) -> bool {
    let exponent = y_abs >> 52;
    // OR-ing in bit 52 accounts for the implicit bit when the fraction is zero.
    // An odd integer's least significant set bit has place value 2^0.
    let trailing_zeros = (y_abs | EXP_MASK).trailing_zeros() as u64;
    exponent + trailing_zeros == 1023 + 52
}

#[cold]
fn pow_exp10_fallback(x: f64) -> f64 {
    f_exp10(x)
}

#[cold]
fn pow_exp2_fallback(x: f64) -> f64 {
    f_exp2(x)
}

#[cold]
#[inline(never)]
fn f_powi(x: f64, y: f64) -> f64 {
    // The caller checks that |y| is an integer in 1..=102.
    let n = unsafe { y.abs().to_int_unchecked::<usize>() };
    // For |x| >= 2^-9, |x|^n >= 2^-918. For smaller |x| in
    // [2^e, 2^(e+1)), the exponent of |x|^n lies in [e*n, (e+1)*n).
    if x.abs() < f64::from_bits(0x3f60_0000_0000_0000) && y.is_sign_positive() {
        let lower_exponent = get_exponent_f64(x) * n as i64;
        if lower_exponent < -968 {
            let negative_result = x.is_sign_negative() && n % 2 != 0;
            if lower_exponent + n as i64 <= -1075 {
                // The magnitude is strictly below half the least subnormal.
                return if negative_result { -0.0 } else { 0.0 };
            }
            return f_powi_underflow(x.abs(), n, negative_result);
        }
    }
    let mut iter_count = n;

    let mut s = DoubleDouble::new(0., x);

    // exponentiation by squaring: O(log(y)) complexity
    let mut acc = if iter_count % 2 != 0 {
        s
    } else {
        DoubleDouble::new(0., 1.)
    };

    while {
        iter_count >>= 1;
        iter_count
    } != 0
    {
        s = DoubleDouble::mult(s, s);
        if iter_count % 2 != 0 {
            acc = DoubleDouble::mult(acc, s);
        }
    }

    let dz = if POWI_HAS_RECIP && y.is_sign_negative() {
        // Below 2^-968, the products forming the accumulator can lose low bits
        // to underflow before we take the reciprocal. A product at or above
        // 2^-968 has operand exponents summing to at least -969, so its exact
        // 106-bit significand fits down to the subnormal grid (2^-1074).
        if acc.hi.abs() < f64::from_bits(0x0370_0000_0000_0000) {
            let sign = if x.is_sign_negative() && is_odd_exponent(y.to_bits() & i64::MAX as u64) {
                -1.0
            } else {
                1.0
            };
            return pow_rational128(x.abs(), y, sign);
        }
        acc.recip()
    } else {
        acc
    };
    let ub = dz.hi + f_fmla(f64::from_bits(0x3c40000000000000), -dz.hi, dz.lo); // 2^-59
    let lb = dz.hi + f_fmla(f64::from_bits(0x3c40000000000000), dz.hi, dz.lo); // 2^-59
    if ub == lb {
        return dz.to_f64();
    }
    f_powi_hard(x, y)
}

#[cold]
#[inline(never)]
fn f_powi_hard(x: f64, y: f64) -> f64 {
    // As in f_powi, |y| is an integer in 1..=102.
    let mut iter_count = unsafe { y.abs().to_int_unchecked::<usize>() };

    let mut s = TripleDouble::new(0., 0., x);

    // exponentiation by squaring: O(log(y)) complexity
    let mut acc = if iter_count % 2 != 0 {
        s
    } else {
        TripleDouble::new(0., 0., 1.)
    };

    while {
        iter_count >>= 1;
        iter_count
    } != 0
    {
        s = TripleDouble::quick_mult(s, s);
        if iter_count % 2 != 0 {
            acc = TripleDouble::quick_mult(acc, s);
        }
    }

    let dz = if POWI_HAS_RECIP && y.is_sign_negative() {
        acc.recip()
    } else {
        acc
    };
    dz.to_f64()
}

#[cold]
#[inline(never)]
fn f_powi_underflow(x: f64, n: usize, negative_result: bool) -> f64 {
    let mut power = DyadicFloat128::new_from_f64(x);
    let mut lower = DyadicFloat128::new_from_f64(1.0);
    let mut remaining = n;
    while remaining != 0 {
        if remaining & 1 != 0 {
            lower = lower * power;
        }
        remaining >>= 1;
        if remaining != 0 {
            power = power * power;
        }
    }

    // Each dyadic multiplication truncates downward by less than 2^-126
    // relative. The expanded product has at most n multiplications, so
    // lower >= x^n * (1 - n*2^-126). For n <= 102 the missing relative
    // amount is below 2^-119 of lower, or 512 units of the mantissa.
    // The binary64 grid has step 2^shift in those same units.
    let shift = (-lower.exponent as i32 - 1074).max(75);
    let spans_midpoint = if shift <= 128 {
        let midpoint = 1u128 << (shift - 1);
        let remainder = lower.mantissa & (midpoint | (midpoint - 1));
        midpoint.wrapping_sub(remainder) <= 512
    } else {
        // The first nonzero result starts just above half a subnormal unit.
        shift == 129 && lower.mantissa > u128::MAX - 512
    };
    let value = if spans_midpoint {
        f_powi_underflow_exact(x, n)
    } else {
        lower
    };
    round_pow_dyadic(value, negative_result)
}

#[cold]
#[inline(never)]
fn f_powi_underflow_exact(x: f64, n: usize) -> DyadicFloat128 {
    // x = m * 2^e with a 53-bit integer significand m. m^n fits in 5406 bits.
    const LIMBS: usize = (53usize * 102).div_ceil(64);
    let (m, exponent) = odd_parts(x);
    let scale = i64::from(exponent) * n as i64;
    let mut value = [0u64; LIMBS];
    value[0] = m;
    let mut len = 1;
    for _ in 1..n {
        let mut carry = 0;
        for limb in &mut value[..len] {
            // The carry is below m, keeping this product below 2^117.
            let product = *limb as u128 * m as u128 + carry;
            *limb = product as u64;
            carry = product >> 64;
        }
        if carry != 0 {
            value[len] = carry as u64;
            len += 1;
        }
    }

    // Keep the leading 128 bits and fold the discarded bits into a sticky
    // bit. The binary64 round bit is at least 74 positions above that bit.
    let bits = 64 * (len - 1) + (64 - value[len - 1].leading_zeros()) as usize;
    let shift = bits.saturating_sub(128);
    let word = shift / 64;
    let offset = shift % 64;
    let pair = (value[word + 1] as u128) << 64 | value[word] as u128;
    let mut top = pair >> offset;
    if offset != 0 {
        top |= (value[word + 2] as u128) << (128 - offset);
    }
    // m is odd, so any nonempty discarded suffix contains a set bit.
    if shift != 0 {
        top |= 1;
    }
    DyadicFloat128::new(DyadicSign::Pos, (scale + shift as i64) as i16, top)
}

#[inline]
fn rsqrt_midpoint_above_root(x: u64, midpoint: u64) -> bool {
    // The midpoint exceeds 1/sqrt(m) exactly when m times its square
    // exceeds 1. Here m and the midpoint are on the 2^-52 and 2^-54
    // grids, so compare the integer product with 2^160. Equality is
    // impossible because the midpoint numerator is odd and exceeds 1.
    let square = u128::from(midpoint) * u128::from(midpoint);
    let x = u128::from(x);
    let product_hi = (square >> 64) * x + ((square as u64 as u128 * x) >> 64);
    product_hi >= 1 << 96
}

#[inline(never)]
fn pow_inverse_sqrt(x: f64) -> f64 {
    let mut bits = x.to_bits();
    let mut e = 0;
    if bits < 1 << 52 {
        bits = (x * f64::from_bits(0x4350_0000_0000_0000)).to_bits();
        e -= 54;
    }
    e += (bits >> 52) as i32 - 1023;
    let significand = bits & ((1 << 52) - 1);
    // Write x = m * 2^(2 * floor(e/2)), with 1 <= m < 4. The final
    // scaling is exact and the normalized result lies in (1/2, 1].
    let m = f64::from_bits(significand | ((1023 + (e & 1)) as u64) << 52);
    // The reciprocal rounding moves sqrt(1/m) by less than 2^-54,
    // since the square-root denominator exceeds 1. Rounding the square
    // root adds at most 2^-54, so r and its two neighbors suffice.
    let r = (1.0 / m).sqrt();
    let r_bits = r.to_bits();
    let scale = f64::from_bits(((1023 - (e >> 1)) as u64) << 52);

    let x_grid = ((1 << 52) | significand) << (e & 1);
    let r_grid = (r_bits - 0.5f64.to_bits()) + (1 << 52);
    let upper = 2 * r_grid + 1;
    let lower = 2 * r_grid - 1;
    let result = if !rsqrt_midpoint_above_root(x_grid, upper) {
        r_bits + 1
    } else if rsqrt_midpoint_above_root(x_grid, lower) {
        r_bits - 1
    } else {
        r_bits
    };
    f64::from_bits(result) * scale
}

/// Raises `x` to the power `y`.
pub fn f_pow(x: f64, y: f64) -> f64 {
    let mut y = y;
    let x_sign = x.is_sign_negative();
    let y_sign = y.is_sign_negative();

    let x_u = x.to_bits();
    let x_a = x_u & 0x7fff_ffff_ffff_ffff;
    let y_a = y.to_bits() & 0x7fff_ffff_ffff_ffff;
    let mut x = x;

    if x.is_nan() || y.is_nan() {
        return f64::NAN;
    }

    let mut s = 1.0;

    // Exceptional exponents.
    if y == 0.0 {
        return 1.0;
    }

    match y_a {
        0x3fe0_0000_0000_0000 => {
            if x == 0.0 || x_u == f64::NEG_INFINITY.to_bits() {
                // pow(-0, 1/2) = +0
                // pow(-inf, 1/2) = +inf
                return if y_sign { 1.0 / (x * x) } else { x * x };
            }
            return if y_sign {
                if x.is_infinite() {
                    return 0.0;
                }
                if x < 0.0 {
                    return f64::NAN;
                }
                #[cfg(any(
                    all(
                        any(target_arch = "x86", target_arch = "x86_64"),
                        target_feature = "fma"
                    ),
                    target_arch = "aarch64"
                ))]
                {
                    let r = x.sqrt() / x;
                    let rx = r * x;
                    let drx = f_fmla(r, x, -rx);
                    let h = f_fmla(r, rx, -1.0) + r * drx;
                    let correction = (-r * 0.5) * h;
                    // The Newton correction and the rounded interval endpoints
                    // are enclosed using 2^-100 * r. See notes/pow_bounds.py.
                    let radius = r * f64::from_bits(0x39b0_0000_0000_0000);
                    let lower = r + (correction - radius);
                    let upper = r + (correction + radius);
                    if lower == upper {
                        return lower;
                    }
                }
                pow_inverse_sqrt(x)
            } else {
                x.sqrt()
            };
        }
        0x3ff0_0000_0000_0000 => {
            return if y_sign { 1.0 / x } else { x };
        }
        0x4000_0000_0000_0000 => {
            if !y_sign {
                // One binary64 multiplication rounds the exact square, including
                // at the subnormal and overflow boundaries.
                return x * x;
            }
            let x_e = get_exponent_f64(x);
            // For |x| >= 2^538, x^-2 <= 2^-1076 and rounds to +0.
            if x_e > 537 {
                return 0.0;
            }
            // not enough precision to make 0.5 ULP for subnormals
            let a_xe = x_e & 0x7fff_ffff_ffff_ffff;
            if a_xe < 70 {
                let x_sqr = DoubleDouble::from_exact_mult(x, x);
                return x_sqr.recip().to_f64();
            }
        }
        _ => {}
    }

    // |y| > |1075 / log2(1 - 2^-53)|.
    if y_a > 0x43d7_4910_d52d_3052 {
        if y_a == f64::INFINITY.to_bits() {
            if x_a == 0x3ff0_0000_0000_0000 {
                // pow(+-1, +-Inf) = 1.0
                return 1.0;
            }

            // pow (|x| < 1, -inf) = +inf
            // pow (|x| < 1, +inf) = 0.0
            // pow (|x| > 1, -inf) = 0.0
            // pow (|x| > 1, +inf) = +inf
            return if (x_a < 1f64.to_bits()) == y_sign {
                f64::INFINITY
            } else {
                0.0
            };
        }
        // x^y will overflow / underflow in double precision.  Set y to a
        // large enough exponent but not too large, so that the computations
        // won't overflow in double precision.
        y = if y_sign {
            f64::from_bits(0xc630000000000000)
        } else {
            f64::from_bits(0x4630000000000000)
        };
    }

    // y is finite and non-zero.

    if x_u == 1f64.to_bits() {
        // pow(1, y) = 1
        return 1.0;
    }

    if x_a == 0 || x_a == f64::INFINITY.to_bits() {
        let magnitude = if (x_a != 0) ^ y_sign {
            f64::INFINITY
        } else {
            0.0
        };
        return if x_sign && is_odd_exponent(y_a) {
            -magnitude
        } else {
            magnitude
        };
    }
    if x == 2.0 {
        return pow_exp2_fallback(y);
    }
    if x == 10.0 {
        return pow_exp10_fallback(y);
    }

    let is_y_integer = is_integer(y);
    // y is integer and in [-102;102] and |x|<2^10

    // The non-FMA reciprocal can overflow in intermediate operations.
    if is_y_integer
        && y_a <= 0x4059800000000000u64
        && x_a <= 0x4090000000000000u64
        && x_a > 0x3cc0_0000_0000_0000
        && (POWI_HAS_RECIP || !y_sign)
    {
        return f_powi(x, y);
    }

    if x_sign {
        if !is_y_integer {
            return f64::NAN;
        }
        x = -x;
        if is_odd_exponent(y_a) {
            s = -1.0;
        }
    }

    // For positive finite x, |ln(x)| <= 1074*ln(2) < 2^10. Thus
    // |y| <= 2^-64 implies |t| < 2^-54 for t = y*ln(x).
    // If t <= 0, exp(t) >= 1+t. Otherwise exp(t)-1 <= t/(1-t) < 2*t.
    // So x^y lies between the rounding midpoints 1-2^-54 and 1+2^-53.
    if y_a <= 0x3bf0_0000_0000_0000 {
        return 1.0;
    }

    let (l, cancel) = pow_log_1(x);
    let r = DoubleDouble::quick_mult_f64(l, y);
    let res = pow_exp_1(r, s);
    static ERR: [u64; 2] = [0x3bf2700000000000, 0x3c55700000000000];
    let res_min = res.hi + f_fmla(f64::from_bits(ERR[cancel as usize]), -res.hi, res.lo);
    let res_max = res.hi + f_fmla(f64::from_bits(ERR[cancel as usize]), res.hi, res.lo);
    if res_min == res_max {
        return res_max;
    }
    pow_rational128(x, y, s)
}

// Nonzero finite input, with its sign handled by the caller.
fn odd_parts(x: f64) -> (u64, i32) {
    let bits = x.to_bits();
    let exp = ((bits >> 52) & 0x7ff) as i32;
    let mantissa = (bits & 0x000f_ffff_ffff_ffff) | (u64::from(exp != 0) << 52);
    let shift = mantissa.trailing_zeros();
    (mantissa >> shift, exp.max(1) - 1075 + shift as i32)
}

// A binary64 value or rounding midpoint has at most 54 significant bits.
#[cold]
#[inline(never)]
fn pow_rounding_boundary(x: f64, y: f64, s: f64) -> Option<f64> {
    let (mut m, mut e) = odd_parts(x);
    let (n, f) = odd_parts(y);
    let (p, e) = if m == 1 {
        if e == 0 {
            return Some(s);
        }
        // x^y is a power of two with an integer exponent iff the shift
        // of n*e discards only zeros. The unsigned product is below 2^64.
        let magnitude = n * u64::from(e.unsigned_abs());
        let exp = if f >= 0 {
            if f > 10 || magnitude > (1075 >> f) {
                1076
            } else {
                magnitude << f
            }
        } else {
            let shift = f.unsigned_abs();
            if shift > magnitude.trailing_zeros() {
                return None;
            }
            magnitude >> shift
        };
        let exp = exp.min(1076) as i32;
        let exp = if (e < 0) ^ y.is_sign_negative() {
            -exp
        } else {
            exp
        };
        (1, exp)
    } else {
        // A negative y leaves an odd denominator. With m >= 3, 3^34 is
        // the last power of the smallest odd base below 2^54. A 53-bit
        // significand also cannot be a perfect 64th power of an odd base.
        if y.is_sign_negative() || !(-5..=5).contains(&f) || n > (34 >> f.max(0)) {
            return None;
        }
        let n = (n << f.max(0)) as u32;
        if f < 0 {
            let roots = f.unsigned_abs();
            let root_degree = 1 << roots;
            // The denominator 2^roots must divide both e and the
            // multiplicity of every odd prime factor of m.
            if e % root_degree != 0 {
                return None;
            }
            e /= root_degree;
            for _ in 0..roots {
                // m fits exactly in binary64, as does an integral root.
                let root = (m as f64).sqrt() as u64;
                if root * root != m {
                    return None;
                }
                m = root;
            }
        }
        let p = m.checked_pow(n)?;
        if p >= 1 << 54 {
            return None;
        }
        (p, e * n as i32)
    };
    let significant_bits = (64 - p.leading_zeros()) as i32;
    let binade = e + significant_bits - 1;
    if binade < -1075 {
        return Some(0.0 * s);
    }
    if binade > 1023 {
        return Some(f64::INFINITY * s);
    }
    // Subnormal values and midpoints lie on the 2^-1075 grid.
    if e < -1075 {
        return None;
    }
    // p is odd. A one-bit reduction is exactly halfway, so increment
    // the truncated value only when it is odd.
    let halve_even = |m: u64| {
        let half = m >> 1;
        half + (half & 1)
    };
    let encoding = if binade >= -1022 {
        let mantissa = if significant_bits == 54 {
            halve_even(p)
        } else {
            p << (53 - significant_bits)
        };
        // Including the implicit bit also carries across a binade when
        // rounding the significand up to 2^53.
        (((binade + 1022) as u64) << 52) + mantissa
    } else if e == -1075 {
        halve_even(p)
    } else {
        p << (e + 1074)
    };
    Some(f64::from_bits(encoding | (s.to_bits() & (1 << 63))))
}

#[cold]
fn pow_rational128(x: f64, y: f64, s: f64) -> f64 {
    let fraction = x.to_bits() & 0x000f_ffff_ffff_ffff;
    // A power of two has at most one bit in its fraction field. An odd
    // factor can only yield a boundary for y between 1/32 and 34.
    if (fraction & fraction.wrapping_sub(1)) == 0 || (1.0 / 32.0..=34.0).contains(&y) {
        if let Some(result) = pow_rounding_boundary(x, y, s) {
            return result;
        }
    }
    let f_y = DyadicFloat128::new_from_f64(y);

    let r = log_dyadic(x) * f_y;
    round_pow_dyadic(exp_dyadic(r), s < 0.0)
}

#[inline(never)]
fn round_pow_dyadic(mut value: DyadicFloat128, negative: bool) -> f64 {
    value.sign = if negative {
        DyadicSign::Neg
    } else {
        DyadicSign::Pos
    };
    value.fast_as_f64()
}

/// Pow for given value for const context.
/// This is simplified version just to make a good approximation on const context.
pub const fn pow(d: f64, n: f64) -> f64 {
    let value = d.abs();

    let r = n * log(value);
    let c = exp(r);
    if n == 0. {
        return 1.;
    }
    if d < 0.0 {
        let y = n as i32;
        if y % 2 == 0 { c } else { -c }
    } else {
        c
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn powf_test() {
        assert!(
            (pow(2f64, 3f64) - 8f64).abs() < 1e-9,
            "Invalid result {}",
            pow(2f64, 3f64)
        );
        assert!(
            (pow(0.5f64, 2f64) - 0.25f64).abs() < 1e-9,
            "Invalid result {}",
            pow(0.5f64, 2f64)
        );
    }

    #[test]
    fn f_pow_test() {
        assert_eq!(f_pow(
             0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000135264499699371,
            -0.5,
        ), 27189929701044785000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000.);
        assert_eq!(f_pow(
            0.000000000000000000000000000000000000000000000000000021798599361155193,
            -2.,
        ),2104470396771397700000000000000000000000000000000000000000000000000000000000000000000000000000000000000000.);
        assert_eq!(
            f_pow(-25192281723900620000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000.,
            -2.),
            0.
        );
        assert_eq!(
            f_pow(0.000000000000000000000000000000000000000000000000000021799650661798696,
                  -2.),
            2104267423084451500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000.
        );
        assert_eq!(
            f_pow(0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000014916691520383755,
            -2.),
            44942267764413600000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000.
        );
        assert_eq!(
            f_pow(
                0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000135264499699371,
                -0.5,
            ),
            27189929701044785000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000.
        );
        assert_eq!(f_pow(1., f64::INFINITY), 1.);
        assert_eq!(f_pow(2., f64::INFINITY), f64::INFINITY);
        assert_eq!(f_pow(f64::INFINITY, -0.5), 0.);
        assert_eq!(
            f_pow(f64::from_bits(0x3fef_ffff_ffff_fffe), -0.5).to_bits(),
            0x3ff0_0000_0000_0001
        );
        assert_eq!(
            f_pow(f64::from_bits(0x000f_ffff_ffff_ffff), -0.5).to_bits(),
            0x5fe0_0000_0000_0001
        );
        assert_eq!(f_pow(f64::MAX, -0.5).to_bits(), 0x1ff0_0000_0000_0000);
        let x = f64::from_bits(0x5ff0_0000_0000_0000); // 2^512
        assert_eq!(f_pow(x, -2.0).to_bits(), 1 << 50);
        let x = f64::from_bits(0x6180_0000_0000_0000); // 2^537
        assert_eq!(f_pow(x, -2.0).to_bits(), 1);
        assert_eq!(f_pow(-0.0, 1.5).to_bits(), 0.0f64.to_bits());
        assert_eq!(f_pow(-0.0, 3.0).to_bits(), (-0.0f64).to_bits());
        assert_eq!(f_pow(f64::NEG_INFINITY, 1.5), f64::INFINITY);
        assert_eq!(
            f_pow(f64::NEG_INFINITY, f64::from_bits(0x43e0_0000_0000_0000)),
            f64::INFINITY
        );
        assert!(
            (f_pow(2f64, 3f64) - 8f64).abs() < 1e-9,
            "Invalid result {}",
            f_pow(2f64, 3f64)
        );
        assert!(
            (f_pow(0.5f64, 2f64) - 0.25f64).abs() < 1e-9,
            "Invalid result {}",
            f_pow(0.5f64, 2f64)
        );
        assert_eq!(f_pow(2.1f64, 2.7f64), 7.412967494768546);
        assert_eq!(f_pow(27., 1. / 3.), 3.);
        assert_eq!(f_pow(3072.0, 34.0).to_bits(), 0x588d_9fe7_7988_1944);
        assert_eq!(f_pow(25.0, 11.5), 11_920_928_955_078_124.0);
        let x = f64::from_bits((1023 - 215) << 52);
        assert_eq!(f_pow(x, 5.0).to_bits(), 0);
        let below_one = f64::from_bits(1.0f64.to_bits() - 1);
        assert_eq!(f_pow(f64::from_bits(1), below_one).to_bits(), 1);
    }

    #[test]
    fn powi_test() {
        assert_eq!(f_pow(f64::from_bits(0x3cc0_0000_0000_0000), 102.), 0.0);
        assert_eq!(f_pow(3., 3.), 27.);
        assert_eq!(f_pow(3., -3.), 1. / 27.);
        assert_eq!(f_pow(3., 102.), 4.638397686588102e48);
        let x = f64::from_bits(0x3f40_0000_0000_0000); // 2^-11
        assert_eq!(f_pow(x, -102.), f64::INFINITY);
        assert_eq!(f_pow(-x, -101.), f64::NEG_INFINITY);
        let x = f64::from_bits(0x3ce4_e87c_573d_94cf);
        assert_eq!(f_pow(x, 22.).to_bits(), 22);

        let half = f64::from_bits(0x3d40_0000_0000_0000); // (2^-43)^25 = 2^-1075
        let above_half = f64::from_bits(half.to_bits() + 1);
        assert_eq!(f_pow(half, 25.).to_bits(), 0);
        assert_eq!(f_pow(-half, 25.).to_bits(), (-0.0f64).to_bits());
        assert_eq!(f_pow(above_half, 25.).to_bits(), 1);
        assert_eq!(f_pow(-above_half, 25.).to_bits(), (1 << 63) | 1);
        assert_eq!(f_pow(0.000000000000011074474670636028, -22.), 10589880229528372000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000.);
    }
}
