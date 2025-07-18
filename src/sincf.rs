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
use crate::common::{f_fmla, f_fmlaf};
use crate::cosf::{sincos_reduce_big, sincos_reduce0, sincos_reduce1};
use crate::double_double::DoubleDouble;
use crate::sinf::SINCOSF_SIN_TABLE;

#[inline(never)]
#[cold]
fn as_sincf_big(x: f32) -> f32 {
    const B: [u64; 4] = [
        0x3f93bd3cc9be45dc,
        0xbf103c1f081b0833,
        0x3e755d3c6fc9ac1f,
        0xbdce1d3ff281b40d,
    ];
    const A: [u64; 4] = [
        0x3fc921fb54442d17,
        0xbf54abbce6256a39,
        0x3ec466bc5a518c16,
        0xbe232bdc61074ff6,
    ];
    let t = x.to_bits();
    let ax = t.wrapping_shl(1);
    if ax >= 0xffu32 << 24 {
        // nan or +-inf
        if ax.wrapping_shl(8) != 0 {
            return x + x;
        }; // nan
        return f32::NAN; // to raise FE_INVALID
    }
    let (z, ia) = sincos_reduce_big(t);
    let z2 = z * z;
    let z4 = z2 * z2;

    let w0 = f_fmla(z2, f64::from_bits(A[1]), f64::from_bits(A[0]));
    let w1 = f_fmla(z2, f64::from_bits(A[3]), f64::from_bits(A[2]));

    let aa = f_fmla(z4, w1, w0);

    let q0 = f_fmla(z2, f64::from_bits(B[1]), f64::from_bits(B[0]));
    let q1 = f_fmla(z2, f64::from_bits(B[3]), f64::from_bits(B[2]));

    let bb = f_fmla(z4, q1, q0);

    let s0 = f64::from_bits(SINCOSF_SIN_TABLE[(ia & 31) as usize]);
    let c0 = f64::from_bits(SINCOSF_SIN_TABLE[((ia.wrapping_add(8)) & 31) as usize]);

    let f0 = f_fmla(-bb, z * s0, aa * c0);
    let r = f_fmla(z, f0, s0);
    let dd = DoubleDouble::from_exact_div(r, x as f64);
    dd.to_f64() as f32
}

/// Computes sinc(x)
///
/// Max found ULP 0.5
#[inline]
pub fn f_sincf(x: f32) -> f32 {
    let t = x.to_bits();
    let ax = t.wrapping_shl(1);
    let ia;
    let z0 = x;
    let z: f64;

    if !x.is_finite() {
        return f32::NAN;
    }

    if x.abs() == 0. {
        return 1.;
    }

    #[allow(clippy::manual_range_contains)]
    if ax > 0x99000000u32 || ax < 0x73000000u32 {
        // |x| > 6.71089e+07 or |x| < 0.000244141
        if ax < 0x73000000u32 {
            // |x| < 0.000244141
            if ax < 0x66000000u32 {
                // |x| < 2.98023e-08
                if ax == 0u32 {
                    return x;
                }
                let sin_x = f_fmlaf(-x, x.abs(), x) as f64;
                let dd = DoubleDouble::from_exact_div(sin_x, x as f64);
                return dd.to_f64() as f32;
            }
            let x = x as f64;
            let sin_x = (-f64::from_bits(0x3fc5555560000000) * x) * (x * x) + x;
            let dd = DoubleDouble::from_exact_div(sin_x, x);
            return dd.to_f64() as f32;
        }
        return as_sincf_big(x);
    }

    const B: [u64; 4] = [
        0x3f93bd3cc9be45dc,
        0xbf103c1f081b0833,
        0x3e755d3c6fc9ac1f,
        0xbdce1d3ff281b40d,
    ];
    const A: [u64; 4] = [
        0x3fc921fb54442d17,
        0xbf54abbce6256a39,
        0x3ec466bc5a518c16,
        0xbe232bdc61074ff6,
    ];

    if ax < 0x822d97c8u32 {
        (z, ia) = sincos_reduce0(z0);
    } else {
        (z, ia) = sincos_reduce1(z0);
    }
    let z2 = z * z;
    let z4 = z2 * z2;

    let w0 = f_fmla(z2, f64::from_bits(A[1]), f64::from_bits(A[0]));
    let w1 = f_fmla(z2, f64::from_bits(A[3]), f64::from_bits(A[2]));

    let aa = f_fmla(z4, w1, w0);

    let q0 = f_fmla(z2, f64::from_bits(B[1]), f64::from_bits(B[0]));
    let q1 = f_fmla(z2, f64::from_bits(B[3]), f64::from_bits(B[2]));

    let bb = f_fmla(z4, q1, q0);

    let s0 = f64::from_bits(SINCOSF_SIN_TABLE[(ia & 31) as usize]);
    let c0 = f64::from_bits(SINCOSF_SIN_TABLE[((ia.wrapping_add(8)) & 31) as usize]);

    let f0 = f_fmla(aa, z * c0, s0);
    let r = f_fmla(-bb, z2 * s0, f0);
    let dd = DoubleDouble::from_exact_div(r, x as f64);
    dd.to_f64() as f32
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_f_sincf() {
        assert_eq!(f_sincf(0.0), 1.0);
        assert_eq!(f_sincf(0.2), 0.99334663);
    }
}
