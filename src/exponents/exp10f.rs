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

pub(crate) struct ExpBReduc {
    pub(crate) hi: f64,
    pub(crate) lo: f64,
}

const MID_BITS: u32 = 5;
const MID_MASK: usize = (1 << MID_BITS) - 1;
const LOG2_B: f64 = f64::from_bits(0x400a934f0979a371) * (1 << MID_BITS) as f64;
const M_LOGB_2_HI: f64 = f64::from_bits(0xbfd34413509f8000) / (1 << MID_BITS) as f64;
const M_LOGB_2_LO: f64 = f64::from_bits(0x3d380433b83b532a) / (1 << MID_BITS) as f64;
const EXP_2_MID: [u64; 32] = [
    0x3ff0000000000000,
    0x3ff059b0d3158574,
    0x3ff0b5586cf9890f,
    0x3ff11301d0125b51,
    0x3ff172b83c7d517b,
    0x3ff1d4873168b9aa,
    0x3ff2387a6e756238,
    0x3ff29e9df51fdee1,
    0x3ff306fe0a31b715,
    0x3ff371a7373aa9cb,
    0x3ff3dea64c123422,
    0x3ff44e086061892d,
    0x3ff4bfdad5362a27,
    0x3ff5342b569d4f82,
    0x3ff5ab07dd485429,
    0x3ff6247eb03a5585,
    0x3ff6a09e667f3bcd,
    0x3ff71f75e8ec5f74,
    0x3ff7a11473eb0187,
    0x3ff82589994cce13,
    0x3ff8ace5422aa0db,
    0x3ff93737b0cdc5e5,
    0x3ff9c49182a3f090,
    0x3ffa5503b23e255d,
    0x3ffae89f995ad3ad,
    0x3ffb7f76f2fb5e47,
    0x3ffc199bdd85529c,
    0x3ffcb720dcef9069,
    0x3ffd5818dcfba487,
    0x3ffdfc97337b9b5f,
    0x3ffea4afa2a490da,
    0x3fff50765b6e4540,
];

pub(crate) const EXP10F_COEFFS: [u64; 5] = [
    0x40026bb1bbb55515,
    0x40053524c73bd3ea,
    0x4000470591dff149,
    0x3ff2bd7c0a9fbc4d,
    0x3fe1429e74a98f43,
];

/// Range reduction function equivalent to exp_b_range_reduc
#[inline]
pub(crate) fn exp_b_range_reduc(x: f32) -> ExpBReduc {
    let xd = x as f64;

    // kd = round(log2(b) * x)
    let kd = (LOG2_B * xd).round();
    let k = kd as i32;

    // hi = floor(kd / 2^MID_BITS)
    let exp_hi = (k.wrapping_shr(MID_BITS) as u64).wrapping_shl(52); // 52 = fraction bits in f64

    // mh = 2^hi * 2^mid
    let mid_index = (k as usize) & MID_MASK;
    let mh_bits = EXP_2_MID[mid_index].wrapping_add(exp_hi);
    let mh = f64::from_bits(mh_bits);

    // dx = x - (hi + mid) * log(2)
    let z0 = f_fmla(kd, M_LOGB_2_HI, xd);
    let dx = f_fmla(kd, M_LOGB_2_LO, z0);

    ExpBReduc { lo: dx, hi: mh }
}

/// Computes exp10 with FMA
///
/// Max found ULP 0.49999508
#[inline]
pub fn f_exp10f(x: f32) -> f32 {
    let x_u = x.to_bits();
    let x_abs = x_u & 0x7fffffff;

    // When |x| >= log10(2^128), or x is nan
    if x_abs >= 0x421a209b {
        // When x < log10(2^-150) or nan
        if x_u > 0xc2349e35 {
            // exp(-Inf) = 0
            if x.is_infinite() {
                return 0.0;
            }
            // exp(nan) = nan
            if x.is_nan() {
                return x;
            }
            return 0.0;
        }
        // x >= log10(2^128) or nan
        if x > 0. && (x_u >= 0x421a209b) {
            // x is +inf or nan
            return x + f32::INFINITY;
        }
    }

    if x_abs <= 0x3b9a209b {
        if x_u == 0xb25e5bd9 {
            // x = -1.2943e-08
            return 1.;
        }
        // |x| < 2^-25
        // 10^x ~ 1 + log(10) * x
        if x_abs <= 0x32800000 {
            return f_fmlaf(x, f32::from_bits(0x40135da2), 1.0);
        }
    }

    // Range reduction: 10^x = 2^(mid + hi) * 10^lo
    //   rr = (2^(mid + hi), lo)
    let rr = exp_b_range_reduc(x);

    // The low part is approximated by a degree-5 minimax polynomial.
    // 10^lo ~ 1 + COEFFS[0] * lo + ... + COEFFS[4] * lo^5
    let lo2 = rr.lo * rr.lo;
    // c0 = 1 + COEFFS[0] * lo
    let c0 = f_fmla(rr.lo, f64::from_bits(EXP10F_COEFFS[0]), 1.0);
    // c1 = COEFFS[1] + COEFFS[2] * lo
    let c1 = f_fmla(
        rr.lo,
        f64::from_bits(EXP10F_COEFFS[2]),
        f64::from_bits(EXP10F_COEFFS[1]),
    );
    // c2 = COEFFS[3] + COEFFS[4] * lo
    let c2 = f_fmla(
        rr.lo,
        f64::from_bits(EXP10F_COEFFS[4]),
        f64::from_bits(EXP10F_COEFFS[3]),
    );
    // p = c1 + c2 * lo^2
    //   = COEFFS[1] + COEFFS[2] * lo + COEFFS[3] * lo^2 + COEFFS[4] * lo^3
    let p = f_fmla(lo2, c2, c1);
    // 10^lo ~ c0 + p * lo^2
    // 10^x = 2^(mid + hi) * 10^lo
    //      ~ mh * (c0 + p * lo^2)
    //      = (mh * c0) + p * (mh * lo^2)
    f_fmla(p, lo2 * rr.hi, c0 * rr.hi) as f32
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_exp10f() {
        assert_eq!(super::f_exp10f(1.), 10.0);
        assert_eq!(super::f_exp10f(2.), 100.0);
        assert_eq!(super::f_exp10f(3.), 1000.0);
    }
}
