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
use crate::dekker::Dekker;
use crate::exp::{EXP_REDUCE_T0, EXP_REDUCE_T1, to_denormal};
use crate::exp2::ldexp;

/// Computes exp10
///
/// Max found ULP 0.5003.
#[inline]
pub fn f_exp10(x: f64) -> f64 {
    let mut ix = x.to_bits();
    let aix = ix & 0x7fff_ffff_ffff_ffff;
    if aix > 0x40734413509f79feu64 {
        // |x| > 0x40734413509f79fe
        if aix > 0x7ff0000000000000u64 {
            return x + x;
        } // nan
        if aix == 0x7ff0000000000000u64 {
            return if (ix >> 63) != 0 { 0.0 } else { x };
        }
        if (ix >> 63) == 0 {
            return f64::from_bits(0x7fe0000000000000) * 2.0; // x > 308.255
        }
        if aix > 0x407439b746e36b52u64 {
            // x < -323.607
            return f64::from_bits(0x0018000000000000) * f64::from_bits(0x3c80000000000000);
        }
    }

    // check x integer to avoid a spurious inexact exception
    if ix.wrapping_shl(16) == 0 && (aix >> 48) <= 0x4036 {
        let kx = x.round_ties_even();
        if kx == x {
            let k = kx as i64;
            if k >= 0 {
                let mut r = 1.0;
                for _ in 0..k {
                    r *= 10.0;
                }
                return r;
            }
        }
    }
    /* avoid spurious underflow: for |x| <= 2.41082e-17
    exp10(x) rounds to 1 to nearest */
    if aix <= 0x3c7bcb7b1526e50eu64 {
        return 1.0 + x; // |x| <= 2.41082e-17
    }
    let t = (f64::from_bits(0x40ca934f0979a371) * x).round_ties_even();
    let jt = t as i64;
    let i1 = jt & 0x3f;
    let i0 = (jt >> 6) & 0x3f;
    let ie = jt >> 12;
    let t00 = EXP_REDUCE_T0[i0 as usize];
    let t01 = EXP_REDUCE_T1[i1 as usize];
    let t0 = Dekker::new(f64::from_bits(t00.0), f64::from_bits(t00.1));
    let t1 = Dekker::new(f64::from_bits(t01.0), f64::from_bits(t01.1));
    let mut tz = Dekker::quick_mult(t0, t1);
    const L0: f64 = f64::from_bits(0x3f13441350800000);
    const L1: f64 = f64::from_bits(0x3d1f79fef311f12b);
    let dx = f_fmla(-L1, t, f_fmla(-L0, t, x));
    let dx2 = dx * dx;

    const CH: [u64; 4] = [
        0x40026bb1bbb55516,
        0x40053524c73cea69,
        0x4000470591fd74e1,
        0x3ff2bd760a1f32a5,
    ];

    let p0 = f_fmla(dx, f64::from_bits(CH[1]), f64::from_bits(CH[0]));
    let p1 = f_fmla(dx, f64::from_bits(CH[3]), f64::from_bits(CH[2]));

    let p = f_fmla(dx2, p1, p0);

    let mut fh = tz.hi;
    let fx = tz.hi * dx;
    let mut fl = f_fmla(fx, p, tz.lo);
    if ix < 0xc0733a7146f72a42u64 {
        // x > -307.653
        fh = ldexp(fh + fl, ie as u64);
    } else {
        // x <= -307.653: exp10(x) < 2^-1022
        ix = 1u64.wrapping_sub(ie as u64).wrapping_shl(52);
        tz = Dekker::from_exact_add(f64::from_bits(ix), fh);
        fl += tz.lo;
        fh = to_denormal(fh + fl);
    }
    fh
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exp10f() {
        assert_eq!(f_exp10(1.), 10.0);
        assert_eq!(f_exp10(2.), 100.0);
        assert_eq!(f_exp10(3.), 1000.0);
        assert_eq!(f_exp10(4.), 10000.0);
        assert_eq!(f_exp10(5.), 100000.0);
        assert_eq!(f_exp10(6.), 1000000.0);
        assert_eq!(f_exp10(7.), 10000000.0);
    }
}
