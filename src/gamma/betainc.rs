/*
 * // Copyright (c) Radzivon Bartoshyk 9/2025. All rights reserved.
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
use crate::bessel::i0_exp;
use crate::common::f_fmla;
use crate::double_double::DoubleDouble;
use crate::gamma::lnbeta::lnbeta_core;
use crate::logs::{fast_log_dd, log1p_fast_dd};
use crate::{f_exp, f_log, f_pow};

/// Regularized incomplete beta
pub fn f_betainc_reg(a: f64, b: f64, x: f64) -> f64 {
    let aa = a.to_bits();
    let ab = b.to_bits();
    let ax = x.to_bits();

    if aa >= 0x7ffu64 << 52
        || aa == 0
        || ab >= 0x7ffu64 << 52
        || ab == 0
        || ax == 0
        || ax >= 0x3ff0000000000000
    {
        if (aa >> 63) != 0 || (ab >> 63) != 0 || (ax >> 63) != 0 {
            // |a| < 0 or |b| < 0
            return f64::NAN;
        }
        if ax >= 0x3ff0000000000000 {
            // |x| > 1
            if ax == 0x3ff0000000000000 {
                // x == 1
                return 1.;
            }
            return f64::NAN;
        }
        if ax.wrapping_shl(1) == 0 {
            // |x| == 0
            return 0.;
        }
        if aa.wrapping_shl(1) == 0 {
            // |a| == 0
            return 1.0;
        }
        if ab.wrapping_shl(1) == 0 {
            // |b| == 0
            return 0.;
        }
        if a.is_infinite() {
            // |a| == inf
            return 0.;
        }
        if b.is_infinite() {
            // |b| == inf
            return 1.;
        }
        return a + f64::NAN; // nan
    }

    if ab == 0x3ff0000000000000 {
        // b == 1
        return f_pow(x, a);
    }

    /*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
    /*Use the fact that beta is symmetrical.*/
    let mut return_inverse = false;
    let mut dx = DoubleDouble::new(0., x);
    let mut a = a;
    let mut b = b;
    if x > (a + 1.0) / (a + b + 2.0) {
        std::mem::swap(&mut a, &mut b);
        dx = DoubleDouble::from_full_exact_sub(1.0, x);
        return_inverse = true;
    }
    /*Find the first part before the continued fraction.*/
    let da = a;
    let db = b;
    let ln_beta_ab = lnbeta_core(a, b);
    let log_dx = fast_log_dd(dx);
    let mut log1p_dx = log1p_fast_dd(-dx.hi);
    log1p_dx.lo += -dx.lo / dx.hi;
    let z1 = DoubleDouble::mul_f64_add(log1p_dx, db, -ln_beta_ab);
    let w0 = DoubleDouble::mul_f64_add(log_dx, da, z1);

    /*Use Lentz's algorithm to evaluate the continued fraction.*/
    let mut f = 1.0;
    let mut c = 1.0;
    let mut d = 0.0;

    const TINY: f64 = 1.0e-31;
    const STOP: f64 = f64::EPSILON;

    for i in 0..200 {
        let m = i / 2;
        let numerator: f64 = if i == 0 {
            1.0 /*First numerator is 1.0.*/
        } else if i % 2 == 0 {
            let m = m as f64;
            let c0 = f_fmla(2.0, m, da);
            (m * (db - m) * dx.hi) / ((c0 - 1.0) * c0) /*Even term.*/
        } else {
            let m = m as f64;
            let c0 = f_fmla(2.0, m, da);
            -((da + m) * (da + db + m) * dx.hi) / (c0 * (c0 + 1.)) /*Odd term.*/
        };

        /*Do an iteration of Lentz's algorithm.*/
        d = f_fmla(numerator, d, 1.0);
        if d.abs() < TINY {
            d = TINY;
        }
        d = 1.0 / d;

        c = 1.0 + numerator / c;
        if c.abs() < TINY {
            c = TINY;
        }

        let cd = c * d;
        f *= cd;

        /*Check for stop.*/
        if (1.0 - cd).abs() < STOP {
            break;
        }
    }

    let r = DoubleDouble::from_full_exact_sub(f, 1.0);
    let exponent = w0.to_f64();
    const LOG_MIN_NORMAL: f64 = -708.3964185322641;
    let exp = if exponent < LOG_MIN_NORMAL {
        DoubleDouble::default()
    } else {
        i0_exp(exponent)
    };
    // Combine in log space when the exponential or its quotient would be
    // subnormal. Keep the residual of the log in the final multiplication.
    if exp.hi <= da * f64::MIN_POSITIVE {
        let log_multiplier = f_log(r.to_f64()) - f_log(da);
        let log_tail = DoubleDouble::quick_dd_add(w0, DoubleDouble::new(0., log_multiplier));
        let tail = f_exp(log_tail.hi);
        let tail = if tail.is_finite() {
            f64::mul_add(tail, log_tail.lo, tail)
        } else {
            tail
        };
        return if return_inverse { 1.0 - tail } else { tail };
    }

    let front = DoubleDouble::div_dd_f64(exp, da);
    if return_inverse {
        DoubleDouble::mul_add_f64(-front, r, 1.).to_f64()
    } else {
        DoubleDouble::quick_mult(front, r).to_f64()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_betainc() {
        assert_eq!(f_betainc_reg(0.5, 2.5, 0.5), 0.9244131815783875);
        assert_eq!(f_betainc_reg(2.5, 1.0, 0.5), 0.1767766952966368811);
        assert_eq!(f_betainc_reg(0.5, 0., 1.), 1.);
        assert_eq!(f_betainc_reg(5., 1.4324, 0.1312), 8.872581630413704e-5);
        assert_eq!(f_betainc_reg(7., 42., 0.4324), 0.9999954480481231);
        assert_eq!(f_betainc_reg(5., 2., 1.), 1.);
        assert_eq!(f_betainc_reg(5., 2., 0.), 0.);
        assert_eq!(f_betainc_reg(5., 2., 0.5), 0.10937500000000006);
        assert!(f_betainc_reg(5., 2., -1.).is_nan());
        assert!(f_betainc_reg(5., 2., 1.1).is_nan());
        assert!(f_betainc_reg(5., 2., f64::INFINITY).is_nan());
        assert!(f_betainc_reg(5., 2., f64::NEG_INFINITY).is_nan());
        assert!(f_betainc_reg(5., 2., f64::NAN).is_nan());
        assert!(f_betainc_reg(-5., 2., 0.432).is_nan());
        assert!(f_betainc_reg(5., -2., 0.432).is_nan());
        assert!(f_betainc_reg(5., 2., -0.432).is_nan());
    }

    #[test]
    fn test_betainc_subnormal_tail() {
        assert_eq!(f_betainc_reg(496.7, 52.0, 0.17), 1.1956759935e-314);
        assert_eq!(f_betainc_reg(510.0, 52.0, 0.17), f64::from_bits(1));
        assert_eq!(f_betainc_reg(511.0, 52.0, 0.17), 0.0);
        assert_eq!(f_betainc_reg(52.0, 496.7, 0.83), 1.0);

        // Dividing the normal exponential by a gives a subnormal quotient.
        let value = f_betainc_reg(1e6, 2.0, 0.9992787037398628);
        let expected = 3.095284631839e-311_f64;
        assert!(
            value.to_bits().abs_diff(expected.to_bits()) <= 4,
            "{value:e}"
        );
    }
}
