/*
 * // Copyright (c) Radzivon Bartoshyk 6/2025. All rights reserved.
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
use crate::acospi::PI_OVER_TWO_F128;
use crate::asin::asin_eval;
use crate::asin_eval_dyadic::asin_eval_dyadic;
use crate::common::f_fmla;
use crate::double_double::DoubleDouble;
use crate::dyadic_float::{DyadicFloat128, DyadicSign};

#[inline]
/// Computes acos(x)
///
/// Max found ULP 0.5
pub fn f_acos(x: f64) -> f64 {
    let x_e = (x.to_bits() >> 52) & 0x7ff;
    const E_BIAS: u64 = (1u64 << (11 - 1u64)) - 1u64;

    const PI_OVER_TWO: DoubleDouble = DoubleDouble::new(
        f64::from_bits(0x3c91a62633145c07),
        f64::from_bits(0x3ff921fb54442d18),
    );

    let x_abs = f64::from_bits(x.to_bits() & 0x7fff_ffff_ffff_ffff);

    // |x| < 0.5.
    if x_e < E_BIAS - 1 {
        // |x| < 2^-55.
        if x_e < E_BIAS - 55 {
            // When |x| < 2^-55, acos(x) = pi/2
            return (x_abs + f64::from_bits(0x35f0000000000000)) + PI_OVER_TWO.hi;
        }

        let x_sq = DoubleDouble::from_exact_mult(x, x);
        let err = x_abs * f64::from_bits(0x3cc0000000000000);
        // Polynomial approximation:
        //   p ~ asin(x)/x
        let (p, err) = asin_eval(x_sq, err);
        // asin(x) ~ x * p
        let r0 = DoubleDouble::from_exact_mult(x, p.hi);
        // acos(x) = pi/2 - asin(x)
        //         ~ pi/2 - x * p
        //         = pi/2 - x * (p.hi + p.lo)
        let r_hi = f_fmla(-x, p.hi, PI_OVER_TWO.hi);
        // Use Dekker's 2SUM algorithm to compute the lower part.
        let mut r_lo = ((PI_OVER_TWO.hi - r_hi) - r0.hi) - r0.lo;
        r_lo = f_fmla(-x, p.lo, r_lo + PI_OVER_TWO.lo);

        let r_upper = r_hi + (r_lo + err);
        let r_lower = r_hi + (r_lo - err);

        if r_upper == r_lower {
            return r_upper;
        }

        // Ziv's accuracy test failed, perform 128-bit calculation.

        // Recalculate mod 1/64.
        let idx = (x_sq.hi * f64::from_bits(0x4050000000000000)).round() as usize;

        // Get x^2 - idx/64 exactly.  When FMA is available, double-double
        // multiplication will be correct for all rounding modes. Otherwise, we use
        // Float128 directly.
        let mut x_f128 = DyadicFloat128::new_from_f64(x);

        let u: DyadicFloat128;
        #[cfg(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        ))]
        {
            // u = x^2 - idx/64
            let u_hi = DyadicFloat128::new_from_f64(f_fmla(
                idx as f64,
                f64::from_bits(0xbf90000000000000),
                x_sq.hi,
            ));
            u = u_hi.quick_add(&DyadicFloat128::new_from_f64(x_sq.lo));
        }

        #[cfg(not(any(
            all(
                any(target_arch = "x86", target_arch = "x86_64"),
                target_feature = "fma"
            ),
            all(target_arch = "aarch64", target_feature = "neon")
        )))]
        {
            let x_sq_f128 = x_f128.quick_mul(&x_f128);
            u = x_sq_f128.quick_add(&DyadicFloat128::new_from_f64(
                idx as f64 * f64::from_bits(0xbf90000000000000),
            ));
        }

        let p_f128 = asin_eval_dyadic(u, idx);
        // Flip the sign of x_f128 to perform subtraction.
        x_f128.sign = x_f128.sign.negate();
        let r = PI_OVER_TWO_F128.quick_add(&x_f128.quick_mul(&p_f128));
        return r.fast_as_f64();
    }

    // |x| >= 0.5

    let x_sign = if x.is_sign_negative() { -1.0 } else { 1.0 };

    const PI: DoubleDouble = DoubleDouble::new(
        f64::from_bits(0x3ca1a62633145c07),
        f64::from_bits(0x400921fb54442d18),
    );

    // |x| >= 1
    if x_e >= E_BIAS {
        // x = +-1, asin(x) = +- pi/2
        if x_abs == 1.0 {
            // x = 1, acos(x) = 0,
            // x = -1, acos(x) = pi
            return if x == 1.0 {
                0.0
            } else {
                f_fmla(-x_sign, PI.hi, PI.lo)
            };
        }
        // |x| > 1, return NaN.
        return f64::NAN;
    }

    // When |x| >= 0.5, we perform range reduction as follow:
    //
    // When 0.5 <= x < 1, let:
    //   y = acos(x)
    // We will use the double angle formula:
    //   cos(2y) = 1 - 2 sin^2(y)
    // and the complement angle identity:
    //   x = cos(y) = 1 - 2 sin^2 (y/2)
    // So:
    //   sin(y/2) = sqrt( (1 - x)/2 )
    // And hence:
    //   y/2 = asin( sqrt( (1 - x)/2 ) )
    // Equivalently:
    //   acos(x) = y = 2 * asin( sqrt( (1 - x)/2 ) )
    // Let u = (1 - x)/2, then:
    //   acos(x) = 2 * asin( sqrt(u) )
    // Moreover, since 0.5 <= x < 1:
    //   0 < u <= 1/4, and 0 < sqrt(u) <= 0.5,
    // And hence we can reuse the same polynomial approximation of asin(x) when
    // |x| <= 0.5:
    //   acos(x) ~ 2 * sqrt(u) * P(u).
    //
    // When -1 < x <= -0.5, we reduce to the previous case using the formula:
    //   acos(x) = pi - acos(-x)
    //           = pi - 2 * asin ( sqrt( (1 + x)/2 ) )
    //           ~ pi - 2 * sqrt(u) * P(u),
    // where u = (1 - |x|)/2.

    // u = (1 - |x|)/2
    let u = f_fmla(x_abs, -0.5, 0.5);
    // v_hi + v_lo ~ sqrt(u).
    // Let:
    //   h = u - v_hi^2 = (sqrt(u) - v_hi) * (sqrt(u) + v_hi)
    // Then:
    //   sqrt(u) = v_hi + h / (sqrt(u) + v_hi)
    //            ~ v_hi + h / (2 * v_hi)
    // So we can use:
    //   v_lo = h / (2 * v_hi).
    let v_hi = u.sqrt();

    let h;
    #[cfg(any(
        all(
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "fma"
        ),
        all(target_arch = "aarch64", target_feature = "neon")
    ))]
    {
        h = f_fmla(v_hi, -v_hi, u);
    }
    #[cfg(not(any(
        all(
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "fma"
        ),
        all(target_arch = "aarch64", target_feature = "neon")
    )))]
    {
        let v_hi_sq = DoubleDouble::from_exact_mult(v_hi, v_hi);
        h = (u - v_hi_sq.hi) - v_hi_sq.lo;
    }

    // Scale v_lo and v_hi by 2 from the formula:
    //   vh = v_hi * 2
    //   vl = 2*v_lo = h / v_hi.
    let vh = v_hi * 2.0;
    let vl = h / v_hi;

    // Polynomial approximation:
    //   p ~ asin(sqrt(u))/sqrt(u)
    let err = vh * f64::from_bits(0x3cc0000000000000);

    let (p, err) = asin_eval(DoubleDouble::new(0.0, u), err);

    // Perform computations in double-double arithmetic:
    //   asin(x) = pi/2 - (v_hi + v_lo) * (ASIN_COEFFS[idx][0] + p)
    let r0 = DoubleDouble::quick_mult(DoubleDouble::new(vl, vh), p);

    let r_hi;
    let r_lo;
    if x.is_sign_positive() {
        r_hi = r0.hi;
        r_lo = r0.lo;
    } else {
        let r = DoubleDouble::from_exact_add(PI.hi, -r0.hi);
        r_hi = r.hi;
        r_lo = (PI.lo - r0.lo) + r.lo;
    }

    let r_upper = r_hi + (r_lo + err);
    let r_lower = r_hi + (r_lo - err);

    if r_upper == r_lower {
        return r_upper;
    }

    // Ziv's accuracy test failed, we redo the computations in Float128.
    // Recalculate mod 1/64.
    let idx = (u * f64::from_bits(0x4050000000000000)).round() as usize;

    // After the first step of Newton-Raphson approximating v = sqrt(u), we have
    // that:
    //   sqrt(u) = v_hi + h / (sqrt(u) + v_hi)
    //      v_lo = h / (2 * v_hi)
    // With error:
    //   sqrt(u) - (v_hi + v_lo) = h * ( 1/(sqrt(u) + v_hi) - 1/(2*v_hi) )
    //                           = -h^2 / (2*v * (sqrt(u) + v)^2).
    // Since:
    //   (sqrt(u) + v_hi)^2 ~ (2sqrt(u))^2 = 4u,
    // we can add another correction term to (v_hi + v_lo) that is:
    //   v_ll = -h^2 / (2*v_hi * 4u)
    //        = -v_lo * (h / 4u)
    //        = -vl * (h / 8u),
    // making the errors:
    //   sqrt(u) - (v_hi + v_lo + v_ll) = O(h^3)
    // well beyond 128-bit precision needed.

    // Get the rounding error of vl = 2 * v_lo ~ h / vh
    // Get full product of vh * vl
    let vl_lo;
    #[cfg(any(
        all(
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "fma"
        ),
        all(target_arch = "aarch64", target_feature = "neon")
    ))]
    {
        vl_lo = f_fmla(-v_hi, vl, h) / v_hi;
    }
    #[cfg(not(any(
        all(
            any(target_arch = "x86", target_arch = "x86_64"),
            target_feature = "fma"
        ),
        all(target_arch = "aarch64", target_feature = "neon")
    )))]
    {
        let vh_vl = DoubleDouble::from_exact_mult(v_hi, vl);
        vl_lo = ((h - vh_vl.hi) - vh_vl.lo) / v_hi;
    }
    let t = h * (-0.25) / u;
    let vll = f_fmla(vl, t, vl_lo);
    let m_v_p = DyadicFloat128::new_from_f64(vl) + DyadicFloat128::new_from_f64(vll);
    let mut m_v = DyadicFloat128::new_from_f64(vh) + m_v_p;
    m_v.sign = if x.is_sign_negative() {
        DyadicSign::Neg
    } else {
        DyadicSign::Pos
    };

    // Perform computations in Float128:
    //   acos(x) = (v_hi + v_lo + vll) * P(u)         , when 0.5 <= x < 1,
    //           = pi - (v_hi + v_lo + vll) * P(u)    , when -1 < x <= -0.5.
    let y_f128 =
        DyadicFloat128::new_from_f64(f_fmla(idx as f64, f64::from_bits(0xbf90000000000000), u));

    let p_f128 = asin_eval_dyadic(y_f128, idx);
    let mut r_f128 = m_v * p_f128;

    if x.is_sign_negative() {
        const PI_F128: DyadicFloat128 = DyadicFloat128 {
            sign: DyadicSign::Pos,
            exponent: -126,
            mantissa: 0xc90fdaa2_2168c234_c4c6628b_80dc1cd1_u128,
        };
        r_f128 = PI_F128 + r_f128;
    }

    r_f128.fast_as_f64()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn f_acos_test() {
        assert_eq!(f_acos(0.7), 0.7953988301841436);
        assert_eq!(f_acos(-0.1), 1.6709637479564565);
        assert_eq!(f_acos(-0.4), 1.9823131728623846);
    }
}
