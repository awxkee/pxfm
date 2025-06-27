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
use crate::common::{f_fmla, f_fmlaf};

static SIN_PI: [u64; 128] = [
    0x0000000000000000,
    0x3fa91f65f10dd814,
    0x3fb917a6bc29b42c,
    0x3fc2c8106e8e613a,
    0x3fc8f8b83c69a60b,
    0x3fcf19f97b215f1b,
    0x3fd294062ed59f06,
    0x3fd58f9a75ab1fdd,
    0x3fd87de2a6aea963,
    0x3fdb5d1009e15cc0,
    0x3fde2b5d3806f63b,
    0x3fe073879922ffee,
    0x3fe1c73b39ae68c8,
    0x3fe30ff7fce17035,
    0x3fe44cf325091dd6,
    0x3fe57d69348ceca0,
    0x3fe6a09e667f3bcd,
    0x3fe7b5df226aafaf,
    0x3fe8bc806b151741,
    0x3fe9b3e047f38741,
    0x3fea9b66290ea1a3,
    0x3feb728345196e3e,
    0x3fec38b2f180bdb1,
    0x3feced7af43cc773,
    0x3fed906bcf328d46,
    0x3fee212104f686e5,
    0x3fee9f4156c62dda,
    0x3fef0a7efb9230d7,
    0x3fef6297cff75cb0,
    0x3fefa7557f08a517,
    0x3fefd88da3d12526,
    0x3feff621e3796d7e,
    0x3ff0000000000000,
    0x3feff621e3796d7e,
    0x3fefd88da3d12526,
    0x3fefa7557f08a517,
    0x3fef6297cff75cb0,
    0x3fef0a7efb9230d7,
    0x3fee9f4156c62dda,
    0x3fee212104f686e5,
    0x3fed906bcf328d46,
    0x3feced7af43cc773,
    0x3fec38b2f180bdb1,
    0x3feb728345196e3e,
    0x3fea9b66290ea1a3,
    0x3fe9b3e047f38741,
    0x3fe8bc806b151741,
    0x3fe7b5df226aafaf,
    0x3fe6a09e667f3bcd,
    0x3fe57d69348ceca0,
    0x3fe44cf325091dd6,
    0x3fe30ff7fce17035,
    0x3fe1c73b39ae68c8,
    0x3fe073879922ffee,
    0x3fde2b5d3806f63b,
    0x3fdb5d1009e15cc0,
    0x3fd87de2a6aea963,
    0x3fd58f9a75ab1fdd,
    0x3fd294062ed59f06,
    0x3fcf19f97b215f1b,
    0x3fc8f8b83c69a60b,
    0x3fc2c8106e8e613a,
    0x3fb917a6bc29b42c,
    0x3fa91f65f10dd814,
    0x0000000000000000,
    0xbfa91f65f10dd814,
    0xbfb917a6bc29b42c,
    0xbfc2c8106e8e613a,
    0xbfc8f8b83c69a60b,
    0xbfcf19f97b215f1b,
    0xbfd294062ed59f06,
    0xbfd58f9a75ab1fdd,
    0xbfd87de2a6aea963,
    0xbfdb5d1009e15cc0,
    0xbfde2b5d3806f63b,
    0xbfe073879922ffee,
    0xbfe1c73b39ae68c8,
    0xbfe30ff7fce17035,
    0xbfe44cf325091dd6,
    0xbfe57d69348ceca0,
    0xbfe6a09e667f3bcd,
    0xbfe7b5df226aafaf,
    0xbfe8bc806b151741,
    0xbfe9b3e047f38741,
    0xbfea9b66290ea1a3,
    0xbfeb728345196e3e,
    0xbfec38b2f180bdb1,
    0xbfeced7af43cc773,
    0xbfed906bcf328d46,
    0xbfee212104f686e5,
    0xbfee9f4156c62dda,
    0xbfef0a7efb9230d7,
    0xbfef6297cff75cb0,
    0xbfefa7557f08a517,
    0xbfefd88da3d12526,
    0xbfeff621e3796d7e,
    0xbff0000000000000,
    0xbfeff621e3796d7e,
    0xbfefd88da3d12526,
    0xbfefa7557f08a517,
    0xbfef6297cff75cb0,
    0xbfef0a7efb9230d7,
    0xbfee9f4156c62dda,
    0xbfee212104f686e5,
    0xbfed906bcf328d46,
    0xbfeced7af43cc773,
    0xbfec38b2f180bdb1,
    0xbfeb728345196e3e,
    0xbfea9b66290ea1a3,
    0xbfe9b3e047f38741,
    0xbfe8bc806b151741,
    0xbfe7b5df226aafaf,
    0xbfe6a09e667f3bcd,
    0xbfe57d69348ceca0,
    0xbfe44cf325091dd6,
    0xbfe30ff7fce17035,
    0xbfe1c73b39ae68c8,
    0xbfe073879922ffee,
    0xbfde2b5d3806f63b,
    0xbfdb5d1009e15cc0,
    0xbfd87de2a6aea963,
    0xbfd58f9a75ab1fdd,
    0xbfd294062ed59f06,
    0xbfcf19f97b215f1b,
    0xbfc8f8b83c69a60b,
    0xbfc2c8106e8e613a,
    0xbfb917a6bc29b42c,
    0xbfa91f65f10dd814,
];

const SIN_COEFFS: [u64; 3] = [0x3da921fb54442d0f, 0xb8f4abbce6102b94, 0x3424669fa3c58463];
const COS_COEFFS: [u64; 3] = [0xbb53bd3cc9be45cf, 0x36903c1f08088742, 0xb1b55d1e5eff55a5];

/// Computes sin(PI*x)
///
/// Max ULP 0.5
#[inline]
pub fn f_sinpif(x: f32) -> f32 {
    let ix = x.to_bits();
    let e: i32 = ((ix >> 23) & 0xff) as i32;
    if e == 0xff {
        if (ix << 9) == 0 {
            return f32::NAN;
        }
        return x + x; // nan
    }
    let mut m: i32 = ((ix & 0x007f_ffff) | (1u32 << 23)) as i32;
    let mut sgn: i32 = ix as i32;
    sgn >>= 31;
    m = (m ^ sgn) - sgn;
    let s = 143i32.wrapping_sub(e);
    if s < 0 {
        // |x| >= 131072
        if s < -6 {
            // |x| >= 8.38861e+06 {
            return f32::copysign(0.0, x);
        }
        let mut iq = (m as u32).wrapping_shl((-s - 1) as u32) as i32;
        iq &= 127;
        if iq == 0 || iq == 64 {
            return f32::copysign(0.0, x);
        }
        return f64::from_bits(SIN_PI[iq as usize]) as f32;
    } else if s > 30 {
        // |x| < 6.10352e-05
        let z = x as f64;
        let z2 = z * z;
        let zw0 = f_fmla(
            f64::from_bits(0xc014abbce625be53),
            z2,
            f64::from_bits(0x400921fb54442d18),
        );
        return (z * zw0) as f32;
    }
    let si = 25i32.wrapping_sub(s);
    if si >= 0 && ((m as u32).wrapping_shl(si as u32)) == 0 {
        return f32::copysign(0.0, x);
    }

    let k = (m as u32).wrapping_shl((31 - s) as u32) as i32;
    let z = k as f64;
    let z2 = z * z;

    let fs0 = f_fmla(
        z2,
        f64::from_bits(SIN_COEFFS[2]),
        f64::from_bits(SIN_COEFFS[1]),
    );

    let fc0 = f_fmla(
        z2,
        f64::from_bits(COS_COEFFS[2]),
        f64::from_bits(COS_COEFFS[1]),
    );

    let fs = f_fmla(z2, fs0, f64::from_bits(SIN_COEFFS[0]));
    let fc = f_fmla(z2, fc0, f64::from_bits(COS_COEFFS[0]));
    let mut iq = (m >> s) as u32;
    iq = (iq.wrapping_add(1)) >> 1;
    let is = iq & 127;
    let ic = (iq + 32) & 127;
    let ts = f64::from_bits(SIN_PI[is as usize]);
    let tc = f64::from_bits(SIN_PI[ic as usize]);

    let r0 = f_fmla(ts * z2, fc, ts);
    let r = f_fmla(tc * z, fs, r0);
    r as f32
}

/// Computes cos(PI*x)
///
/// Max ULP 0.5
#[inline]
pub fn f_cospif(x: f32) -> f32 {
    let ix = x.to_bits();
    let e: i32 = ((ix >> 23) & 0xff) as i32;
    if e == 0xff {
        if (ix.wrapping_shl(9)) == 0 {
            return f32::NAN;
        }
        return x + x; // nan
    }
    let m: i32 = ((ix & 0x007f_ffff) | (1u32 << 23)) as i32;
    let s: i32 = 143i32.wrapping_sub(e);
    let p: i32 = e.wrapping_sub(112);
    if p < 0
    // |x| < 2^-15
    {
        let ax: u32 = ix & 0x7fffffff;
        // Warning: -4.9348 * x underflows for |x| < 2.38205e-39
        if ax >= 0x19f030u32 {
            return f_fmlaf(f32::from_bits(0xc09de9e6) * x, x, 1.0);
        } else {
            // |x| < 2.38205e-39
            return f_fmlaf(-x, x, 1.0);
        }
    }
    if p > 31 {
        if p > 63 {
            return 1.0;
        }
        let iq: i32 = (m as u32).wrapping_shl((p - 32) as u32) as i32;
        return f64::from_bits(SIN_PI[((iq.wrapping_add(32)) & 127) as usize]) as f32;
    }
    let k: i32 = (m as u32).wrapping_shl(p as u32) as i32;
    if k == 0 {
        let iq = m >> (32u32.wrapping_sub(p as u32));
        return f64::from_bits(SIN_PI[((iq.wrapping_add(32)) & 127) as usize]) as f32;
    }
    let z = k as f64;
    let z2 = z * z;
    let fs0 = f_fmla(
        z2,
        f64::from_bits(SIN_COEFFS[2]),
        f64::from_bits(SIN_COEFFS[1]),
    );

    let fc0 = f_fmla(
        z2,
        f64::from_bits(COS_COEFFS[2]),
        f64::from_bits(COS_COEFFS[1]),
    );

    let fs = f_fmla(z2, fs0, f64::from_bits(SIN_COEFFS[0]));
    let fc = f_fmla(z2, fc0, f64::from_bits(COS_COEFFS[0]));

    let mut iq: u32 = (m >> s) as u32;
    iq = (iq.wrapping_add(1)) >> 1;
    let is: u32 = iq & 127;
    let ic: u32 = (iq + 32) & 127;
    let ts = f64::from_bits(SIN_PI[ic as usize]);
    let tc = f64::from_bits(SIN_PI[is as usize]);

    let r0 = f_fmla(ts * z2, fc, ts);
    let r = f_fmla(-(tc * z), fs, r0);

    r as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_f_sinpif() {
        assert_eq!(f_sinpif(115.30706), -0.82185423);
        assert!(f_sinpif(f32::INFINITY).is_nan());
    }

    #[test]
    fn test_f_cospif() {
        assert_eq!(f_sinpif(115.30706), -0.5696978);
        assert!(f_cospif(f32::INFINITY).is_nan());
    }
}
