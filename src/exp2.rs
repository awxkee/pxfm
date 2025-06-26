/*
 * // Copyright (c) Radzivon Bartoshyk 4/2025. All rights reserved.
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
use crate::dekker::Dekker;
use crate::exp::{EXP_REDUCE_T0, EXP_REDUCE_T1, to_denormal};

#[inline]
pub(crate) fn ldexp(d: f64, i: u64) -> f64 {
    let b = d.to_bits();
    f64::from_bits(b.wrapping_add(i.wrapping_shl(52)))
}

/// Computes exp2
///
/// Max found ULP 0.5015
#[inline]
pub fn f_exp2(x: f64) -> f64 {
    let mut ix = x.to_bits();
    let ax = ix.wrapping_shl(1);
    if ax == 0 {
        return 1.0;
    }
    if ax >= 0x8120000000000000u64 {
        // |x| >= 1024
        if ax > 0xffe0000000000000u64 {
            return x + x; // nan
        }
        if ax == 0xffe0000000000000u64 {
            return if (ix >> 63) != 0 { 0.0 } else { x };
        }
        // +/-inf
        if (ix >> 63) != 0 {
            // x <= -1024
            if ix >= 0xc090cc0000000000u64 {
                // x <= -1075
                const Z: f64 = f64::from_bits(0x0010000000000000);
                return Z * Z;
            }
        } else {
            // x >= 1024
            return f64::from_bits(0x7fe0000000000000) * x;
        }
    }

    // for |x| <= 0x1.71547652b82fep-54, 2^x rounds to 1 to nearest
    // this avoids a spurious underflow in z*z below
    if ax <= 0x792e2a8eca5705fcu64 {
        return 1.0 + f64::copysign(f64::from_bits(0x3c90000000000000), x);
    }

    let m = ix.wrapping_shl(12);
    let ex = (ax >> 53).wrapping_sub(0x3ff);
    let frac = ex >> 63 | m << (ex & 63);
    let sx = 4096.0 * x;
    let fx = sx.round_ties_even();
    let z = sx - fx;
    let z2 = z * z;
    let k = fx as i64;
    let i1 = k & 0x3f;
    let i0 = (k >> 6) & 0x3f;
    let ie = k >> 12;
    let t00 = EXP_REDUCE_T0[i0 as usize];
    let t01 = EXP_REDUCE_T1[i1 as usize];
    let t0 = Dekker::new(f64::from_bits(t00.0), f64::from_bits(t00.1));
    let t1 = Dekker::new(f64::from_bits(t01.0), f64::from_bits(t01.1));
    let ti0 = Dekker::quick_mult(t0, t1);
    const C: [u64; 4] = [
        0x3f262e42fefa39ef,
        0x3e4ebfbdff82c58f,
        0x3d6c6b08d73b3e01,
        0x3c83b2ab6fdda001,
    ];
    let tz = ti0.hi * z;
    let mut fh = ti0.hi;

    let p0 = f_fmla(z, f64::from_bits(C[1]), f64::from_bits(C[0]));
    let p1 = f_fmla(z, f64::from_bits(C[3]), f64::from_bits(C[2]));
    let p2 = f_fmla(z2, p1, p0);

    let mut fl = f_fmla(tz, p2, ti0.lo);

    const EPS: f64 = f64::from_bits(0x3c0833beace2b6fe);

    if ix <= 0xc08ff00000000000u64 {
        // x >= -1022
        if frac != 0 {
            fh += fl - EPS;
        }
        fh = ldexp(fh, ie as u64);
    } else {
        // subnormal case
        ix = 1u64.wrapping_sub(ie as u64).wrapping_shl(52);
        let rs = Dekker::from_exact_add(f64::from_bits(ix), fh);
        fl += rs.lo;
        fh = rs.hi;
        if frac != 0 {
            fh += fl - EPS;
        }
        // when 2^x is exact, no underflow should be raised
        fh = to_denormal(fh);
    }
    fh
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exp2d() {
        assert_eq!(f_exp2(2.0), 4.0);
        assert_eq!(f_exp2(3.0), 8.0);
        assert_eq!(f_exp2(4.0), 16.0);
        assert!((f_exp2(0.35f64) - 0.35f64.exp2()).abs() < 1e-8);
        assert!((f_exp2(-0.6f64) - (-0.6f64).exp2()).abs() < 1e-8);
    }
}
