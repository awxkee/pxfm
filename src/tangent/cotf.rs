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
use crate::tangent::tanf::{tan_reduce_big, tan_reduce1};

/// Computes cotangent
///
/// Max found ULP 0.4999999
#[inline]
pub fn f_cotf(x: f32) -> f32 {
    let t = x.to_bits();
    let e = (t.wrapping_shr(23)) & 0xff;
    let i;
    let z;
    if e < 127 + 28 {
        // |x| < 2^28
        if e < 115 {
            // |x| < 2^-13
            if e < 102 {
                // |x| < 2^-26
                return (1. / x as f64) as f32;
            }
            let ddx = x as f64;
            let dx = 1. / ddx;
            // taylor order 3
            return f_fmla(ddx, f64::from_bits(0xbfd5555555555555), dx) as f32;
        } else if e < 117 {
            // |x| < 2^-8
            let ddx = x as f64;
            let dx = 1. / ddx;
            let t2 = ddx * ddx;
            // taylor order 1
            let mut z = f_fmla(
                t2,
                f64::from_bits(0xbf61566abc011567),
                f64::from_bits(0xbf96c16c16c16c17),
            );
            z = f_fmla(t2, z, f64::from_bits(0xbfd5555555555555));
            return f_fmla(z, ddx, dx) as f32;
        } else if e < 119 {
            // |x| < 2^-6
            let ddx = x as f64;
            let dx = 1. / ddx;
            let t2 = ddx * ddx;
            // taylor order 9
            let mut z = f_fmla(
                t2,
                f64::from_bits(0xbef66a8f2bf70ebe),
                f64::from_bits(0xbf2bbd779334ef0b),
            );
            z = f_fmla(t2, z, f64::from_bits(0xbf61566abc011567));
            z = f_fmla(t2, z, f64::from_bits(0xbf96c16c16c16c17));
            z = f_fmla(t2, z, f64::from_bits(0xbfd5555555555555));
            return f_fmla(z, ddx, dx) as f32;
        }
        (z, i) = tan_reduce1(x);
    } else if e < 0xff {
        (z, i) = tan_reduce_big(t);
    } else {
        if (t.wrapping_shl(9)) != 0 {
            return x + x;
        } // nan
        return f32::NAN; // inf
    }
    let z2 = z * z;
    let z4 = z2 * z2;
    const CN: [u64; 4] = [
        0x3ff921fb54442d18,
        0xbfdfd226e573289f,
        0x3f9b7a60c8dac9f6,
        0xbf2725beb40f33e5,
    ];
    const CD: [u64; 4] = [
        0x3ff0000000000000,
        0xbff2395347fb829d,
        0x3fc2313660f29c36,
        0xbf69a707ab98d1c1,
    ];
    // cos/sin
    const S: [f64; 2] = [0., 1.];
    let mut n = f_fmla(z2, f64::from_bits(CN[1]), f64::from_bits(CN[0]));
    let n2 = f_fmla(z2, f64::from_bits(CN[3]), f64::from_bits(CN[2]));
    n = f_fmla(z4, n2, n);
    let mut d = f_fmla(z2, f64::from_bits(CD[1]), f64::from_bits(CD[0]));
    let d2 = f_fmla(z2, f64::from_bits(CD[3]), f64::from_bits(CD[2]));
    d = f_fmla(z4, d2, d);
    n *= z;
    let s0 = S[(i & 1) as usize];
    let s1 = S[(1 - (i & 1)) as usize];
    let r1 = f_fmla(n, s0, d * s1) / f_fmla(n, s1, -d * s0);
    r1 as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cotf_test() {
        assert_eq!(f_cotf(0.0010348097), 966.36084);
        assert_eq!(f_cotf(0.0020286469), 492.93872);
        assert_eq!(f_cotf(-0.0020286469), -492.93872);
        assert_eq!(f_cotf(1.0020286469), 0.63923126);
        assert_eq!(f_cotf(-1.0020286469), -0.63923126);
        assert_eq!(f_cotf(0.0), f32::INFINITY);
        assert_eq!(f_cotf(-0.0), f32::NEG_INFINITY);
        assert!(f_cotf(f32::INFINITY).is_nan());
        assert!(f_cotf(f32::NEG_INFINITY).is_nan());
        assert!(f_cotf(f32::NAN).is_nan());
    }
}
