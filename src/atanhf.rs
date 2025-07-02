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
use crate::asinhf::log_eval;
use crate::common::f_fmla;

#[inline(always)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn f_polyeval6(x: f64, a0: f64, a1: f64, a2: f64, a3: f64, a4: f64, a5: f64) -> f64 {
    let t0 = f_fmla(x, a5, a4); // a3 * x + a2
    let t2 = f_fmla(x, t0, a3); // a3 * x + a2
    let t3 = f_fmla(x, t2, a2); // (a3 * x + a2) * x + a1
    let t4 = f_fmla(x, t3, a1); // (a3 * x + a2) * x + a1
    f_fmla(x, t4, a0) // ((a3 * x + a2) * x + a1) * x + a0
}

/// Hyperbolic atan
///
/// Max ULP 0.5
#[inline]
pub fn f_atanhf(x: f32) -> f32 {
    let x_abs = x.to_bits() & 0x7fff_ffff;

    // |x| >= 1.0
    if x_abs >= 0x3F80_0000u32 {
        if x.is_nan() {
            return x;
        }
        // |x| == 1.0
        return if x_abs == 0x3F80_0000u32 {
            if x.is_sign_positive() {
                f32::INFINITY
            } else {
                f32::NEG_INFINITY
            }
        } else {
            f32::NAN
        };
    }

    // |x| < ~0.10
    if x_abs <= 0x3dcc_0000u32 {
        // |x| <= 2^-26
        if x_abs <= 0x3280_0000u32 {
            return if x_abs == 0 {
                x
            } else {
                (x as f64 + f64::from_bits(0x3fd5555555555555) * x as f64 * x as f64 * x as f64)
                    as f32
            };
        }

        let xdbl = x as f64;
        let x2 = xdbl * xdbl;
        // Pure Taylor series.
        let pe = f_polyeval6(
            x2,
            0.0,
            f64::from_bits(0x3fd5555555555555),
            f64::from_bits(0x3fc999999999999a),
            f64::from_bits(0x3fc2492492492492),
            f64::from_bits(0x3fbc71c71c71c71c),
            f64::from_bits(0x3fb745d1745d1746),
        );
        return f_fmla(xdbl, pe, xdbl) as f32;
    }
    let xdbl = x as f64;
    (0.5 * log_eval((xdbl + 1.0) / (xdbl - 1.0))) as f32
}
