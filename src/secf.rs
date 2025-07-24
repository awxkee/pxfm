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
use crate::polyeval::{f_polyeval4, f_polyeval6};
use std::hint::black_box;

fn as_secf_big(x: f32) -> f32 {
    let t = x.to_bits();
    let ax = t.wrapping_shl(1);
    if ax >= 0xffu32 << 24 {
        if ax << 8 != 0 {
            return x + x;
        }

        return f32::NAN;
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
    let (z, ia) = sincos_reduce_big(t);
    let z2 = z * z;
    let z4 = z2 * z2;

    let w0 = f_fmla(z2, f64::from_bits(A[1]), f64::from_bits(A[0]));
    let w1 = f_fmla(z2, f64::from_bits(A[3]), f64::from_bits(A[2]));

    let aa = f_fmla(z4, w1, w0);

    let q0 = f_fmla(z2, f64::from_bits(B[1]), f64::from_bits(B[0]));
    let q1 = f_fmla(z2, f64::from_bits(B[3]), f64::from_bits(B[2]));

    let bb = f_fmla(z4, q1, q0);

    let s0 = f64::from_bits(crate::cosf::SINCOS_F_TABLE[((ia.wrapping_add(8i32)) & 31) as usize]);
    let c0 = f64::from_bits(crate::cosf::SINCOS_F_TABLE[(ia & 31) as usize]);

    let g0 = f_fmla(aa, s0, -bb * (z * c0));

    let r = 1. / f_fmla(z, g0, c0);
    r as f32
}

/// Computes secant
///
/// Max found ULP 0.5
#[inline]
pub fn f_secf(x: f32) -> f32 {
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
    let (z, ia);
    let z0 = x;
    if ax > 0x99000000u32 || ax < 0x73000000 {
        if ax < 0x73000000 {
            if ax < 0x66000000u32 {
                if ax == 0u32 {
                    return 1.0;
                };
                return black_box(1.0) - black_box(f64::from_bits(0x3e60000000000000) as f32);
            }
            return f_fmlaf(f32::from_bits(0x3f000000) * x, x, 1.0);
        }
        return as_secf_big(x);
    }

    if ax < 0x79eb851eu32 {
        // 0.03
        let dx = x as f64;
        let x2 = dx * dx;
        // taylor order 6
        let p = f_polyeval4(
            x2,
            1.,
            f64::from_bits(0x3fe0000000000000),
            f64::from_bits(0x3fcaaaaaaaaaaaab),
            f64::from_bits(0x3fb5b05b05b05b06),
        );
        return p as f32;
    } else if ax < 0x7bae147au32 {
        // 0.105
        let dx = x as f64;
        let x2 = dx * dx;
        // taylor order 10
        let p = f_polyeval6(
            x2,
            1.,
            f64::from_bits(0x3fe0000000000000),
            f64::from_bits(0x3fcaaaaaaaaaaaab),
            f64::from_bits(0x3fb5b05b05b05b06),
            f64::from_bits(0x3fa1965965965966),
            f64::from_bits(0x3f8c834283cd3723),
        );
        return p as f32;
    }

    if ax < 0x82a41896u32 {
        // Exception
        if ax == 0x7f955ffcu32 {
            return f32::from_bits(0xc29d7d8b);
        }
        // 13.128001
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

    let c0 = f64::from_bits(crate::cosf::SINCOS_F_TABLE[(ia & 31) as usize]);
    let s0 = f64::from_bits(crate::cosf::SINCOS_F_TABLE[(ia.wrapping_add(8) & 31) as usize]);

    let n0 = f_fmla(bb, -(z2 * c0), c0);

    let r = 1. / f_fmla(aa, z * s0, n0);
    r as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_f_secf() {
        assert_eq!(f_secf(0.0), 1.0);
        assert_eq!(f_secf(0.5), 1.139494);
        assert_eq!(f_secf(-0.5), 1.139494);
        assert_eq!(f_secf(1.5), 14.136833);
        assert_eq!(f_secf(-1.5), 14.136833);
        assert!(f_secf(f32::INFINITY).is_nan());
        assert!(f_secf(f32::NEG_INFINITY).is_nan());
        assert!(f_secf(f32::NAN).is_nan());
    }
}
