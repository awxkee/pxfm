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
use crate::common::f_fmla;

/// Computes tan(PI*x)
///
/// Max found ULP 0.5
#[inline]
pub fn f_tanpif(x: f32) -> f32 {
    let mut ix = x.to_bits();
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
        return f32::copysign(0.0, x);
    }
    let x4 = 4.0 * x;
    let nx4 = x4.round_ties_even();
    let dx4 = x4 - nx4;
    let ni = x.round_ties_even();
    let zf = x - ni;
    if dx4 == 0.0 {
        // 4*x integer
        let mut k = x4 as i32;
        if (k & 1) != 0 {
            return f32::copysign(1.0, zf);
        } // x = 1/4 mod 1/2
        k &= 6;
        if k == 0 {
            return f32::copysign(0.0, x);
        } // x = 0 mod 2
        if k == 4 {
            return -f32::copysign(0.0, x);
        } // x = 1 mod 2

        if k == 2 {
            return f32::INFINITY;
        } // x = 1/2 mod 2
        // now necessarily k=6
        return f32::NEG_INFINITY; // x = -1/2 mod 2
    }

    ix = zf.to_bits();
    let a: u32 = ix & 0x7fffffff;
    // x=0.287537 is not correctly rounded for RNDZ/RNDD by the code below
    if a == 0x3e933802u32 {
        return f32::copysign(f32::from_bits(0x3fa267dd), zf)
            + f32::copysign(f32::from_bits(0x33000000), zf);
    }
    // x=-0.000115586 is not correctly rounded for RNDU by the code below
    if a == 0x38f26685u32 {
        return f32::copysign(f32::from_bits(0x39be6182), zf)
            + f32::copysign(f32::from_bits(0x2d000000), zf);
    }
    let z = zf as f64;
    let z2 = z * z;

    const CN: [u64; 4] = [
        0x3fe921fb54442d19,
        0xbfd1f458b3e1f8d6,
        0x3f968a34bd0b8f6a,
        0xbf2e4866f7a25f99,
    ];
    const CD: [u64; 4] = [
        0x3ff0000000000000,
        0xbfe4b4b98d2df3a7,
        0x3fb8e9926d2bb901,
        0xbf6a6f77fd847ee0,
    ];
    let z4 = z2 * z2;

    let den0 = f_fmla(z2, f64::from_bits(CD[3]), f64::from_bits(CD[2]));
    let den1 = f_fmla(z2, f64::from_bits(CD[1]), f64::from_bits(CD[0]));

    let num0 = f_fmla(z2, f64::from_bits(CN[1]), f64::from_bits(CN[0]));
    let num1 = f_fmla(z2, f64::from_bits(CN[3]), f64::from_bits(CN[2]));

    let den = f_fmla(z4, den0, den1) * (0.25 - z2);

    let r0 = f_fmla(-z, z2, z);
    let r1 = f_fmla(z4, num1, num0);

    let r = r0 * r1 / den;
    r as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tanpif() {
        assert_eq!(f_tanpif(115.30706), 1.4426143);
        assert!(f_tanpif(f32::INFINITY).is_nan());
    }
}
