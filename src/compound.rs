/*
 * // Copyright (c) Radzivon Bartoshyk 7/2025. All rights reserved.
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
use crate::common::dd_fmla;
use crate::dekker::Dekker;
use crate::dyadic_float::{DyadicFloat128, DyadicSign};
use crate::log1p::{log1p_f64_dd, log1p_f64_dyadic};
use crate::pow::{is_integer, is_odd_integer};
use crate::pow_exec::{exp_dyadic, pow_exp_1};

/// Computes (1+x)^y
///
/// max found ULP 0.5
#[inline]
pub fn f_compound(x: f64, y: f64) -> f64 {
    /* Rules from IEEE 754-2019 for compound (x, n) with n integer:
       (a) compound (x, 0) is 1 for x >= -1 or quiet NaN
       (b) compound (-1, n) is +Inf and signals the divideByZero exception for n < 0
       (c) compound (-1, n) is +0 for n > 0
       (d) compound (+/-0, n) is 1
       (e) compound (+Inf, n) is +Inf for n > 0
       (f) compound (+Inf, n) is +0 for n < 0
       (g) compound (x, n) is qNaN and signals the invalid exception for x < -1
       (h) compound (qNaN, n) is qNaN for n <> 0.
    */

    let mut y = y;
    let x_sign = x.is_sign_negative();
    let y_sign = y.is_sign_negative();

    let x_abs = x.to_bits() & 0x7fff_ffff_ffff_ffff;
    let y_abs = y.to_bits() & 0x7fff_ffff_ffff_ffff;

    const MANTISSA_MASK: u64 = (1u64 << 52) - 1;

    let y_mant = y.to_bits() & MANTISSA_MASK;
    let x_u = x.to_bits();
    let x_a = x_abs;
    let y_a = y_abs;

    let mut x = x;

    // If x or y is signaling NaN
    if x.is_nan() || y.is_nan() {
        return f64::NAN;
    }

    let mut s = 1.0;

    // The double precision number that is closest to 1 is (1 - 2^-53), which has
    //   log2(1 - 2^-53) ~ -1.715...p-53.
    // So if |y| > |1075 / log2(1 - 2^-53)|, and x is finite:
    //   |y * log2(x)| = 0 or > 1075.
    // Hence, x^y will either overflow or underflow if x is not zero.
    if y_mant == 0
        || y_a > 0x43d7_4910_d52d_3052
        || x_u == 1f64.to_bits()
        || x_u >= f64::INFINITY.to_bits()
        || x_u < f64::MIN.to_bits()
    {
        // Exceptional exponents.
        if y == 0.0 {
            return 1.0;
        }

        // (h) compound(qNaN, n) is qNaN for n ≠ 0
        if x.is_nan() {
            if y != 0. {
                return x;
            } // propagate qNaN
            return 1.0;
        }

        // (d) compound(±0, n) is 1
        if x == 0.0 {
            return 1.0;
        }

        // (e, f) compound(+Inf, n)
        if x.is_infinite() && x > 0.0 {
            return if y > 0. { x } else { 0.0 };
        }

        // (g) compound(x, n) is qNaN and signals invalid for x < -1
        if x < -1.0 {
            // Optional: raise invalid explicitly
            return f64::NAN;
        }

        // (b, c) compound(-1, n)
        if x == -1.0 {
            return if y < 0. { f64::INFINITY } else { 0.0 };
        }

        match y_a {
            0x3fe0_0000_0000_0000 => {
                // TODO: speed up x^(-1/2) with rsqrt(x) when available.
                if x == 0.0 {
                    return 1.0;
                }
                let z = Dekker::from_full_exact_add(x, 1.0).sqrt();
                return if y_sign {
                    z.recip().to_f64()
                } else {
                    z.to_f64()
                };
            }
            0x3ff0_0000_0000_0000 => {
                return if y_sign {
                    Dekker::from_full_exact_add(x, 1.0).recip().to_f64()
                } else {
                    Dekker::from_full_exact_add(x, 1.0).to_f64()
                };
            }
            0x4000_0000_0000_0000 => {
                let z0 = Dekker::from_full_exact_add(x, 1.0);
                let z = Dekker::quick_mult(z0, z0);
                return if y_sign {
                    z.recip().to_f64()
                } else {
                    f64::copysign(z.to_f64(), x)
                };
            }
            _ => {}
        }

        // |y| > |1075 / log2(1 - 2^-53)|.
        if y_a > 0x43d7_4910_d52d_3052 {
            if y_a >= 0x7ff0_0000_0000_0000 {
                // y is inf or nan
                if y_mant != 0 {
                    // y is NaN
                    // pow(1, NaN) = 1
                    // pow(x, NaN) = NaN
                    return if x_u == 1f64.to_bits() { 1.0 } else { y };
                }

                // Now y is +-Inf
                if f64::from_bits(x_abs).is_nan() {
                    // pow(NaN, +-Inf) = NaN
                    return x;
                }

                if x_a == 0x3ff0_0000_0000_0000 {
                    // pow(+-1, +-Inf) = 1.0
                    return 1.0;
                }

                if x == 0.0 && y_sign {
                    // pow(+-0, -Inf) = +inf and raise FE_DIVBYZERO
                    return f64::INFINITY;
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

        if x == 0.0 {
            let out_is_neg = x_sign && is_odd_integer(y);
            if y_sign {
                // pow(0, negative number) = inf
                return if out_is_neg {
                    f64::NEG_INFINITY
                } else {
                    f64::INFINITY
                };
            }
            // pow(0, positive number) = 0
            return if out_is_neg { -0.0 } else { 0.0 };
        }

        if x_a == f64::INFINITY.to_bits() {
            let out_is_neg = x_sign && is_odd_integer(y);
            if y_sign {
                return if out_is_neg { -0.0 } else { 0.0 };
            }
            return if out_is_neg {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            };
        }

        if x_a > f64::INFINITY.to_bits() {
            // x is NaN.
            // pow (aNaN, 0) is already taken care above.
            return x;
        }

        // x is finite and negative, and y is a finite integer.
        if x_sign {
            if is_integer(y) {
                x = -x;
                if is_odd_integer(y) {
                    // sign = -1.0;
                    static CS: [f64; 2] = [1.0, -1.0];

                    // set sign to 1 for y even, to -1 for y odd
                    let y_parity = if (y.abs()) >= f64::from_bits(0x4340000000000000) {
                        0usize
                    } else {
                        (y as i64 & 0x1) as usize
                    };
                    s = CS[y_parity];
                }
            } else {
                // pow( negative, non-integer ) = NaN
                return f64::NAN;
            }
        }

        // y is finite and non-zero.

        if x_u == 1f64.to_bits() {
            // pow(1, y) = 1
            return 1.0;
        }

        if x == 0.0 {
            let out_is_neg = x_sign && is_odd_integer(y);
            if y_sign {
                // pow(0, negative number) = inf
                return if out_is_neg {
                    f64::NEG_INFINITY
                } else {
                    f64::INFINITY
                };
            }
            // pow(0, positive number) = 0
            return if out_is_neg { -0.0 } else { 0.0 };
        }

        if x_a == f64::INFINITY.to_bits() {
            let out_is_neg = x_sign && is_odd_integer(y);
            if y_sign {
                return if out_is_neg { -0.0 } else { 0.0 };
            }
            return if out_is_neg {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            };
        }

        if x_a > f64::INFINITY.to_bits() {
            // x is NaN.
            // pow (aNaN, 0) is already taken care above.
            return x;
        }
    }

    let ax = x.to_bits() & 0x7fff_ffff_ffff_ffff;
    let ay = y.to_bits() & 0x7fff_ffff_ffff_ffff;

    // evaluate (1+x)^y explicitly for integer y in [-16,16] range and |x|<2^64
    if y.floor() == y && ay <= 0x4030_0000_0000_0000u64 && ax <= 0x43e0_0000_0000_0000u64 {
        if ax <= 0x3cc0_0000_0000_0000 {
            let mut p = Dekker::from_exact_mult(y, x);
            p = Dekker::full_add_f64(p, 1.0);
            return p.to_f64();
        } // does it work for |x|<2^-29 and |y|<=16?
        let z0 = f64::from_bits(ay) as f32;
        let ay0 = z0.to_bits().wrapping_shl(1);
        let ky: i32 = (((ay0 & 0x00ffffff) | 1 << 24) >> (151 - (ay0 >> 24))) as i32;
        let s = 1.0 + x;
        let mut p = 1.;
        let s2 = s * s;
        let s4 = s2 * s2;
        let s8 = s4 * s4;
        let s16 = s8 * s8;
        let sn: [f64; 6] = [1., s, s2, s4, s8, s16];
        p *= sn[(ky & 1) as usize];
        p *= sn[(ky & 2) as usize];
        p *= sn[(((ky >> 2) & 1) * 3) as usize];
        p *= sn[((ky >> 1) & 4) as usize];
        p *= sn[(((ky >> 4) & 1) * 5) as usize];
        return if y.is_sign_negative() {
            Dekker::new(0., p).recip().to_f64()
        } else {
            p
        };
    }

    // approximate log(x)
    let mut l = log1p_f64_dd(x);

    /* We should avoid a spurious underflow/overflow in y*log(x).
    Underflow: for x<>1, the smallest absolute value of log(x) is obtained
    for x=1-2^-53, with |log(x)| ~ 2^-53. Thus to avoid a spurious underflow
    we require |y| >= 2^-969.
    Overflow: the largest absolute value of log(x) is obtained for x=2^-1074,
    with |log(x)| < 745. Thus to avoid a spurious overflow we require
    |y| < 2^1014. */
    let ey = ((y.to_bits() >> 52) & 0x7ff) as i32;
    if ey < 0x36 || ey >= 0x7f5 {
        l.lo = f64::NAN;
        l.hi = f64::NAN;
    }

    let r = Dekker::quick_mult_f64(l, y);
    if r.hi < 300. {
        let res = pow_exp_1(r, s);

        let res_min = res.hi + dd_fmla(f64::from_bits(0x3c55700000000000), -res.hi, res.lo);
        let res_max = res.hi + dd_fmla(f64::from_bits(0x3c55700000000000), res.hi, res.lo);
        if res_min == res_max {
            return res_max;
        }
    }

    compound_accurate(x, y, s)
}

#[cold]
fn compound_accurate(x: f64, y: f64, s: f64) -> f64 {
    /* the idea of returning res_max instead of res_min is due to Laurent
    Théry: it is better in case of underflow since res_max = +0 always. */

    let f_y = DyadicFloat128::new_from_f64(y);

    let r = log1p_f64_dyadic(x) * f_y;
    let mut result = exp_dyadic(r);

    // 2^R.ex <= R < 2^(R.ex+1)

    /* case R < 2^-1075: underflow case */
    if result.exponent + 127 < -1075 {
        return 0.5 * (s * f64::from_bits(0x0000000000000001));
    }

    result.sign = if s == -1.0 {
        DyadicSign::Neg
    } else {
        DyadicSign::Pos
    };

    result.fast_as_f64()
}

#[cfg(test)]
mod tests {
    use crate::f_compound;

    #[test]
    fn test_compound() {
        assert_eq!(f_compound(2., 5.), 243.);
        assert_eq!(f_compound(126.4324324, 126.4324324), 1.4985383310514043e266);
        assert_eq!(f_compound(0.4324324, 126.4324324), 5.40545942023447e19);
        assert!(f_compound(-0.4324324, 126.4324324).is_nan());
        assert_eq!(f_compound(0.0, 0.0), 1.0);
        assert_eq!(f_compound(0.0, -1. / 2.), 1.0);
        assert_eq!(f_compound(-1., -1. / 2.), f64::INFINITY);
        assert_eq!(f_compound(f64::INFINITY, -1. / 2.), 0.0);
        assert_eq!(f_compound(f64::INFINITY, 1. / 2.), f64::INFINITY);
        assert_eq!(f_compound(46.3828125, 46.3828125), 5.248159634773675e77);
    }
}
