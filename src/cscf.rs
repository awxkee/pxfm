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
use crate::common::f_fmla;
use crate::cosf::{sincos_reduce_big, sincos_reduce0, sincos_reduce1};
use crate::polyeval::f_polyeval3;
use crate::sinf::SINCOSF_SIN_TABLE;

fn as_cscf_big(x: f32) -> f32 {
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
        return f32::NAN;
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
    let r = 1. / f_fmla(z, f0, s0);
    r as f32
}

/// Cosecant
///
/// Max found ULP 0.4999996
#[inline]
pub fn f_cscf(x: f32) -> f32 {
    let t = x.to_bits();
    let ax = t.wrapping_shl(1);
    let ia;
    let z0 = x;
    let z: f64;
    #[allow(clippy::manual_range_contains)]
    if ax > 0x99000000u32 || ax < 0x73000000u32 {
        // |x| > 6.71089e+07 or |x| < 0.000244141
        if ax < 0x73000000u32 {
            // |x| < 0.000244141
            if ax < 0x66000000u32 {
                // |x| < 2.98023e-08
                if ax == 0u32 {
                    return if x.is_sign_negative() {
                        f32::NEG_INFINITY
                    } else {
                        f32::INFINITY
                    };
                }
                let dx = x as f64;
                return (1. / dx) as f32;
            }
            let dx = x as f64;
            let c_term = 1. / dx;
            return f_fmla(dx, f64::from_bits(0x3fc5555555555555), c_term) as f32;
        }
        return as_cscf_big(x);
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
        if ax < 0x79eb851eu32 {
            // 0.03
            let dx = x as f64;
            let c_term = 1. / dx;
            let x2 = dx * dx;
            let p = f_polyeval3(
                x2,
                f64::from_bits(0x3fc5555555555555),
                f64::from_bits(0x3f93e93e93e93e94),
                f64::from_bits(0x3f60b2463814bc5f),
            );
            return f_fmla(dx, p, c_term) as f32;
        }

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
    let r = 1. / f_fmla(-bb, z2 * s0, f0);
    r as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f_cscf_test() {
        assert_eq!(f_cscf(0.04915107), 20.353632);
        assert_eq!(f_cscf(0.5), 2.0858297);
        assert_eq!(f_cscf(0.07), 14.297387);
        assert_eq!(f_cscf(3.6171106e-5), 27646.375);
        assert_eq!(f_cscf(-5.535772e-10), -1806432800.0);
        assert_eq!(f_cscf(0.0), f32::INFINITY);
        assert_eq!(f_cscf(-0.0), f32::NEG_INFINITY);
        assert_eq!(f_cscf(-1.0854926e-19), -9.2124077e18);
    }
}
