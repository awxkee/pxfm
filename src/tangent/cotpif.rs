/*
 * // Copyright (c) Radzivon Bartoshyk 8/2025. All rights reserved.
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

use crate::common::f_fmla;
use crate::sin_cosf::{ArgumentReducerPi, sincospif_eval, sincospif_eval_argument};

/// Computes 1/tan(PI*x)
///
/// Max found ULP 0.5
#[inline]
pub fn f_cotpif(x: f32) -> f32 {
    let ix = x.to_bits();
    let e = ix & (0xff << 23);
    if e > (150 << 23) {
        // |x| > 2^23
        if e == (0xff << 23) {
            // x = nan or inf
            if (ix.wrapping_shl(9)) == 0 {
                // x = inf
                return f32::NAN;
            }
            return x + x; // x = nan
        }
        return f32::INFINITY;
    }
    let argument_reduction = ArgumentReducerPi { x: x as f64 };

    let (y, k) = argument_reduction.reduce();

    if y == 0.0 {
        let km = (k.abs() & 31) as i32; // k mod 32

        match km {
            0 => return f32::copysign(f32::INFINITY, x), // cotpi(n) = ∞
            16 => return 0.0f32.copysign(x),             // cotpi(n+0.5) = 0
            8 => return f32::copysign(1.0, x),           // cotpi(n+0.25) = 1
            24 => return -f32::copysign(1.0, x),         // cotpi(n+0.75) = -1
            _ => {}
        }
    }

    let rs = sincospif_eval_argument(y, k);
    // tan(x) = sin(x) / cos(x)
    //        = (sin_y * cos_k + cos_y * sin_k) / (cos_y * cos_k - sin_y * sin_k)
    let v_cos = f_fmla(rs.sin_y, -rs.sin_k, f_fmla(rs.cosm1_y, rs.cos_k, rs.cos_k));
    let v_sin = f_fmla(rs.sin_y, rs.cos_k, f_fmla(rs.cosm1_y, rs.sin_k, rs.sin_k));

    (v_cos / v_sin) as f32
}

#[inline]
pub(crate) fn cotpif_core(x: f64) -> f64 {
    let rs = sincospif_eval(x);
    // tan(x) = sin(x) / cos(x)
    //        = (sin_y * cos_k + cos_y * sin_k) / (cos_y * cos_k - sin_y * sin_k)
    let v_cos = f_fmla(rs.sin_y, -rs.sin_k, f_fmla(rs.cosm1_y, rs.cos_k, rs.cos_k));
    let v_sin = f_fmla(rs.sin_y, rs.cos_k, f_fmla(rs.cosm1_y, rs.sin_k, rs.sin_k));

    v_cos / v_sin
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cotpif() {
        assert_eq!(f_cotpif(10775313000000000000000000000000.), f32::INFINITY);
        assert_eq!(f_cotpif(5.5625), -0.19891237);
        assert_eq!(f_cotpif(-29.75), 1.0);
        assert_eq!(f_cotpif(-21.5625), 0.19891237);
        assert_eq!(f_cotpif(-15.611655), 0.3659073);
        assert_eq!(f_cotpif(115.30706), 0.693186);
        assert_eq!(f_cotpif(0.), f32::INFINITY);
        assert!(f_cotpif(f32::INFINITY).is_nan());
        assert!(f_cotpif(f32::NAN).is_nan());
    }
}
