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
use crate::acospi::INV_PI_DD;
use crate::atan2::{ATAN_I, atan_eval};
use crate::common::f_fmla;
use crate::dekker::Dekker;

/// Computes atan(x)/PI
///
/// Max found ULP 0.5
#[inline]
pub fn f_atan2pi(y: f64, x: f64) -> f64 {
    const IS_NEG: [f64; 2] = [1.0, -1.0];
    const ZERO: Dekker = Dekker::new(0.0, 0.0);
    const MZERO: Dekker = Dekker::new(-0.0, -0.0);
    const PI: Dekker = Dekker::new(
        f64::from_bits(0x3ca1a62633145c07),
        f64::from_bits(0x400921fb54442d18),
    );
    const MPI: Dekker = Dekker::new(
        f64::from_bits(0xbca1a62633145c07),
        f64::from_bits(0xc00921fb54442d18),
    );
    const PI_OVER_2: Dekker = Dekker::new(
        f64::from_bits(0x3c91a62633145c07),
        f64::from_bits(0x3ff921fb54442d18),
    );
    const MPI_OVER_2: Dekker = Dekker::new(
        f64::from_bits(0xbc91a62633145c07),
        f64::from_bits(0xbff921fb54442d18),
    );
    const PI_OVER_4: Dekker = Dekker::new(
        f64::from_bits(0x3c81a62633145c07),
        f64::from_bits(0x3fe921fb54442d18),
    );
    const THREE_PI_OVER_4: Dekker = Dekker::new(
        f64::from_bits(0x3c9a79394c9e8a0a),
        f64::from_bits(0x4002d97c7f3321d2),
    );

    // Adjustment for constant term:
    //   CONST_ADJ[x_sign][y_sign][recip]
    const CONST_ADJ: [[[Dekker; 2]; 2]; 2] = [
        [[ZERO, MPI_OVER_2], [MZERO, MPI_OVER_2]],
        [[MPI, PI_OVER_2], [MPI, PI_OVER_2]],
    ];

    let x_sign = x.is_sign_negative() as usize;
    let y_sign = y.is_sign_negative() as usize;
    let x_bits = x.to_bits() & 0x7fff_ffff_ffff_ffff;
    let y_bits = y.to_bits() & 0x7fff_ffff_ffff_ffff;
    let x_abs = x_bits;
    let y_abs = y_bits;
    let recip = x_abs < y_abs;
    let mut min_abs = if recip { x_abs } else { y_abs };
    let mut max_abs = if !recip { x_abs } else { y_abs };
    let mut min_exp = min_abs.wrapping_shr(52);
    let mut max_exp = max_abs.wrapping_shr(52);

    let mut num = f64::from_bits(min_abs);
    let mut den = f64::from_bits(max_abs);

    // Check for exceptional cases, whether inputs are 0, inf, nan, or close to
    // overflow, or close to underflow.
    if max_exp > 0x7ffu64 - 128u64 || min_exp < 128u64 {
        if x.is_nan() || y.is_nan() {
            return f64::NAN;
        }
        let x_except = if x == 0.0 {
            0
        } else if x.is_infinite() {
            2
        } else {
            1
        };
        let y_except = if y == 0.0 {
            0
        } else if y.is_infinite() {
            2
        } else {
            1
        };

        // Exceptional cases:
        //   EXCEPT[y_except][x_except][x_is_neg]
        // with x_except & y_except:
        //   0: zero
        //   1: finite, non-zero
        //   2: infinity
        const EXCEPTS: [[[Dekker; 2]; 3]; 3] = [
            [[ZERO, PI], [ZERO, PI], [ZERO, PI]],
            [[PI_OVER_2, PI_OVER_2], [ZERO, ZERO], [ZERO, PI]],
            [
                [PI_OVER_2, PI_OVER_2],
                [PI_OVER_2, PI_OVER_2],
                [PI_OVER_4, THREE_PI_OVER_4],
            ],
        ];

        if (x_except != 1) || (y_except != 1) {
            let mut r = EXCEPTS[y_except][x_except][x_sign];
            r = Dekker::quick_mult(r, INV_PI_DD);
            return f_fmla(IS_NEG[y_sign], r.hi, IS_NEG[y_sign] * r.lo);
        }
        let scale_up = min_exp < 128u64;
        let scale_down = max_exp > 0x7ffu64 - 128u64;
        // At least one input is denormal, multiply both numerator and denominator
        // by some large enough power of 2 to normalize denormal inputs.
        if scale_up {
            num *= f64::from_bits(0x43f0000000000000);
            if !scale_down {
                den *= f64::from_bits(0x43f0000000000000)
            }
        } else if scale_down {
            den *= f64::from_bits(0x3bf0000000000000);
            if !scale_up {
                num *= f64::from_bits(0x3bf0000000000000);
            }
        }

        min_abs = num.to_bits();
        max_abs = den.to_bits();
        min_exp = min_abs.wrapping_shr(52);
        max_exp = max_abs.wrapping_shr(52);
    }
    let final_sign = IS_NEG[if (x_sign != y_sign) != recip { 1 } else { 0 }];
    let const_term = CONST_ADJ[x_sign][y_sign][if recip { 1 } else { 0 }];
    let exp_diff = max_exp - min_exp;
    // We have the following bound for normalized n and d:
    //   2^(-exp_diff - 1) < n/d < 2^(-exp_diff + 1).
    if exp_diff > 54 {
        let r = f_fmla(
            final_sign,
            const_term.hi,
            final_sign * (const_term.lo + num / den),
        );
        let p = Dekker::f64_mult(r, INV_PI_DD);
        return p.to_f64();
    }

    let mut k = (64.0 * num / den).round();
    let idx = k as u64;
    // k = idx / 64
    k *= f64::from_bits(0x3f90000000000000);

    // Range reduction:
    // atan(n/d) - atan(k/64) = atan((n/d - k/64) / (1 + (n/d) * (k/64)))
    //                        = atan((n - d * k/64)) / (d + n * k/64))
    let num_k = Dekker::from_exact_mult(num, k);
    let den_k = Dekker::from_exact_mult(den, k);

    // num_dd = n - d * k
    let num_dd = Dekker::from_exact_add(num - den_k.hi, -den_k.lo);
    // den_dd = d + n * k
    let mut den_dd = Dekker::from_exact_add(den, num_k.hi);
    den_dd.lo += num_k.lo;

    // q = (n - d * k) / (d + n * k)
    let q = Dekker::div(num_dd, den_dd);
    // p ~ atan(q)
    let p = atan_eval(q);

    let vl = ATAN_I[idx as usize];
    let vlo = Dekker::new(f64::from_bits(vl.0), f64::from_bits(vl.1));
    let mut r = Dekker::add(const_term, Dekker::add(vlo, p));
    r = Dekker::quick_mult(r, INV_PI_DD);
    r.hi *= final_sign;
    r.lo *= final_sign;

    r.hi + r.lo
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atan2pi() {
        assert_eq!(f_atan2pi(-5., 2.), -0.3788810584091566);
        assert_eq!(f_atan2pi(2., -5.), 0.8788810584091566);
    }
}
