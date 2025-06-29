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
use crate::atan::{ATAN_CIRCLE, ATAN_REDUCE, poly_dd_3};
use crate::common::f_fmla;
use crate::dekker::Dekker;

const ONE_OVER_PIH: f64 = f64::from_bits(0x3fd45f306dc9c883);
const ONE_OVER_PIL: f64 = f64::from_bits(0xbc76b01ec5417056);
const ONE_OVER_3PI: f64 = f64::from_bits(0x3fbb2995e7b7b604); // approximates 1/(3pi)

#[inline]
fn atanpi_small(x: f64) -> f64 {
    if x == 0. {
        return x;
    }
    if x.abs() == f64::from_bits(0x0015cba89af1f855) {
        return if x > 0. {
            f_fmla(
                f64::from_bits(0x9a70000000000000),
                f64::from_bits(0x1a70000000000000),
                f64::from_bits(0x0006f00f7cd3a40b),
            )
        } else {
            f_fmla(
                f64::from_bits(0x1a70000000000000),
                f64::from_bits(0x1a70000000000000),
                f64::from_bits(0x8006f00f7cd3a40b),
            )
        };
    }
    // generic worst case
    let mut v = x.to_bits();
    if (v & 0xfffffffffffff) == 0x59af9a1194efe
    // +/-0x1.59af9a1194efe*2^e
    {
        let e = v >> 52;
        if (e & 0x7ff) > 2 {
            v = ((e - 2) << 52) | 0xb824198b94a89;
            return if x > 0. {
                f_fmla(
                    f64::from_bits(0x9a70000000000000),
                    f64::from_bits(0x1a70000000000000),
                    f64::from_bits(v),
                )
            } else {
                f_fmla(
                    f64::from_bits(0x1a70000000000000),
                    f64::from_bits(0x1a70000000000000),
                    f64::from_bits(v),
                )
            };
        }
    }
    let h = x * ONE_OVER_PIH;
    /* Assuming h = x*ONE_OVER_PIH - e, the correction term is
    e + x * ONE_OVER_PIL, but we need to scale values to avoid underflow. */
    let mut corr = f_fmla(
        x * f64::from_bits(0x4690000000000000),
        ONE_OVER_PIH,
        -h * f64::from_bits(0x4690000000000000),
    );
    corr = f_fmla(x * f64::from_bits(0x4690000000000000), ONE_OVER_PIL, corr);
    // now return h + corr * 2^-106
    let res = f_fmla(corr, f64::from_bits(0x3950000000000000), h);
    res
}

/* Deal with the case where |x| is large:
for x > 0, atanpi(x) = 1/2 - 1/pi * 1/x + 1/(3pi) * 1/x^3 + O(1/x^5)
for x < 0, atanpi(x) = -1/2 - 1/pi * 1/x + 1/(3pi) * 1/x^3 + O(1/x^5).
The next term 1/5*x^5/pi is smaller than 2^-107 * atanpi(x)
when |x| > 0x1.bep20. */
#[inline]
fn atanpi_asympt(x: f64) -> f64 {
    let h = f64::copysign(0.5, x);
    // approximate 1/x as yh + yl
    let yh = 1.0 / x;
    // Newton's iteration for the inverse is y = y + y*(1-x*y)
    let yl = yh * f_fmla(yh, -x, 1.0);
    let mut m = Dekker::mult(Dekker::new(yl, yh), Dekker::new(ONE_OVER_PIH, ONE_OVER_PIL));
    // m + l ~ 1/pi * 1/x
    m.hi = -m.hi;
    m.lo = f_fmla(ONE_OVER_3PI * yh, yh * yh, -m.lo);
    // m + l ~ - 1/pi * 1/x + 1/(3pi) * 1/x^3
    let vh = Dekker::from_exact_add(h, m.hi);
    m.hi = vh.hi;
    m = Dekker::from_exact_add(vh.lo, m.lo);
    if m.hi.abs() == f64::from_bits(0x3c80000000000000) {
        // this is 1/2 ulp(atan(x))
        m.hi = if m.hi * m.lo > 0. {
            f64::copysign(f64::from_bits(0x3c80000000000001), m.hi)
        } else {
            f64::copysign(f64::from_bits(0x3c7fffffffffffff), m.hi)
        };
    }
    h + m.hi
}

#[inline]
fn atanpi_tiny(x: f64) -> f64 {
    let h = x * ONE_OVER_PIH;
    let mut l = f_fmla(x, ONE_OVER_PIH, -h);
    l = f_fmla(x, ONE_OVER_PIL, l);
    l = f_fmla(-ONE_OVER_3PI * x, x * x, l);
    h + l
}

fn as_atan_refine2(x: f64, a: f64) -> f64 {
    if x.abs() > f64::from_bits(0x413be00000000000) {
        return atanpi_asympt(x);
    }
    if x.abs() < f64::from_bits(0x3e4c700000000000) {
        return atanpi_tiny(x);
    }
    const CH: [(u64, u64); 3] = [
        (0xbfd5555555555555, 0xbc75555555555555),
        (0x3fc999999999999a, 0xbc6999999999bcb8),
        (0xbfc2492492492492, 0xbc6249242093c016),
    ];
    const CL: [u64; 4] = [
        0x3fbc71c71c71c71c,
        0xbfb745d1745d1265,
        0x3fb3b13b115bcbc4,
        0xbfb1107c41ad3253,
    ];
    let phi = ((a.abs()) * f64::from_bits(0x40545f306dc9c883) + 256.5).to_bits();
    let i: i64 = ((phi >> (52 - 8)) & 0xff) as i64;
    let (h, hl);
    if i == 128 {
        h = -1.0 / x;
        hl = f_fmla(h, x, 1.) * h;
    } else {
        let ta = f64::copysign(f64::from_bits(ATAN_REDUCE[i as usize].0), x);
        let zta = x * ta;
        let ztal = f_fmla(x, ta, -zta);
        let zmta = x - ta;
        let v = 1. + zta;
        let d = 1. - v;
        let ev = (d + zta) - ((d + v) - 1.) + ztal;
        let r = 1.0 / v;
        let rl = f_fmla(-ev, r, f_fmla(r, -v, 1.0)) * r;
        h = r * zmta;
        hl = f_fmla(rl, zmta, f_fmla(r, zmta, -h));
    }
    let d2 = Dekker::mult(Dekker::new(hl, h), Dekker::new(hl, h));
    let h4 = d2.hi * d2.hi;
    let h3 = Dekker::mult(Dekker::new(hl, h), d2);

    let fl0 = f_fmla(d2.hi, f64::from_bits(CL[1]), f64::from_bits(CL[0]));
    let fl1 = f_fmla(d2.hi, f64::from_bits(CL[3]), f64::from_bits(CL[2]));

    let fl = d2.hi * f_fmla(h4, fl1, fl0);
    let mut f = poly_dd_3(d2, CH, fl);
    f = Dekker::mult(h3, f);
    let (ah, mut al, mut at);
    if i == 0 {
        ah = h;
        al = f.hi;
        at = f.lo;
    } else {
        let mut df = 0.;
        if i < 128 {
            df = f64::copysign(1.0, x) * f64::from_bits(ATAN_REDUCE[i as usize].1);
        }
        let id = f64::copysign(i as f64, x);
        ah = f64::from_bits(0x3f8921fb54442d00) * id;
        al = f64::from_bits(0x3c88469898cc5180) * id;
        at = f64::from_bits(0xb97fc8f8cbb5bf80) * id;
        let v0 = Dekker::add(Dekker::new(at, al), Dekker::new(0., df));
        let v1 = Dekker::add(v0, Dekker::new(hl, h));
        let v2 = Dekker::add(v1, f);
        al = v2.hi;
        at = v2.lo;
    }

    let v2 = Dekker::from_exact_add(ah, al);
    let v1 = Dekker::from_exact_add(v2.lo, at);

    let z0 = Dekker::mult(
        Dekker::new(v1.hi, v2.hi),
        Dekker::new(ONE_OVER_PIL, ONE_OVER_PIH),
    );
    // atanpi_end
    z0.to_f64()
}

/// Computes atan(x) / pi
#[inline]
pub fn f_atanpi(x: f64) -> f64 {
    const CH: [u64; 4] = [
        0x3ff0000000000000,
        0xbfd555555555552b,
        0x3fc9999999069c20,
        0xbfc248d2c8444ac6,
    ];
    let t = x.to_bits();
    let at: u64 = t & 0x7fff_ffff_ffff_ffff;
    let mut i = (at >> 51).wrapping_sub(2030u64);
    if at < 0x3f7b21c475e6362au64 {
        // |x| < 0.006624
        if at < 0x3c90000000000000u64 {
            // |x| < 2^-54
            return atanpi_small(x);
        }
        if x == 0. {
            return x;
        }
        const CH2: [u64; 4] = [
            0xbfd5555555555555,
            0x3fc99999999998c1,
            0xbfc249249176aec0,
            0x3fbc711fd121ae80,
        ];
        let x2 = x * x;
        let x3 = x * x2;
        let x4 = x2 * x2;
        let f = x3
            * ((f64::from_bits(CH2[0]) + x2 * f64::from_bits(CH2[1]))
                + x4 * (f64::from_bits(CH2[2]) + x2 * f64::from_bits(CH2[3])));
        // begin_atanpi
        /* Here x+f approximates atan(x), with absolute error bounded by
        0x4.8p-52*f (see atan.c). After multiplying by 1/pi this error
        will be bounded by 0x1.6fp-52*f. For |x| < 0x1.b21c475e6362ap-8
        we have |f| < 2^-16*|x|, thus the error is bounded by
        0x1.6fp-52*2^-16*|x| < 0x1.6fp-68. */
        // multiply x + f by 1/pi
        let hy = Dekker::quick_mult(Dekker::new(f, x), Dekker::new(ONE_OVER_PIL, ONE_OVER_PIH));
        /* The rounding error in muldd and the approximation error between
        1/pi and ONE_OVER_PIH + ONE_OVER_PIL are covered by the difference
        between 0x4.8p-52*pi and 0x1.6fp-52, which is > 2^-61.8. */
        let mut ub = hy.hi + f_fmla(f64::from_bits(0x3bb6f00000000000), x, hy.lo);
        let lb = hy.hi + f_fmla(f64::from_bits(0xbbb6f00000000000), x, hy.lo);
        if ub == lb {
            return ub;
        }
        // end_atanpi
        ub = (f + f * f64::from_bits(0x3cd2000000000000)) + x; // atanpi_specific, original value in atan.c
        return as_atan_refine2(x, ub);
    }
    // now |x| >= 0x1.b21c475e6362ap-8
    let h;
    let mut a: Dekker;
    if at > 0x4062ded8e34a9035u64 {
        // |x| > 0x1.2ded8e34a9035p+7, atanpi|x| > 0.49789
        if at >= 0x43445f306dc9c883u64 {
            // |x| >= 0x1.45f306dc9c883p+53, atanpi|x| > 0.5 - 0x1p-55
            if at >= (0x7ffu64 << 52) {
                // case Inf or NaN
                if at == 0x7ffu64 << 52 {
                    // Inf
                    return f64::copysign(0.5, x);
                } // atanpi_specific
                return x + x; // NaN
            }
            return f64::copysign(0.5, x) - f64::copysign(f64::from_bits(0x3c70000000000000), x);
        }
        h = -1.0 / x;
        a = Dekker::new(
            f64::copysign(f64::from_bits(0x3c91a62633145c07), x),
            f64::copysign(f64::from_bits(0x3ff921fb54442d18), x),
        );
    } else {
        // 0x1.b21c475e6362ap-8 <= |x| <= 0x1.2ded8e34a9035p+7
        /* we need to deal with |x| = 1 separately since in this case
        h=0 below, and the error is measured in terms of multiple of h */
        if at == 0x3ff0000000000000 {
            // |x| = 1
            return f64::copysign(f64::from_bits(0x3fd0000000000000), x);
        }
        let u: u64 = t & 0x0007ffffffffffff;
        let ut = u >> (51 - 16);
        let ut2 = ut * ut >> 16;
        let vc = ATAN_CIRCLE[i as usize];
        i = (((vc[0] as u64).wrapping_shl(16)) + ut * (vc[1] as u64) - ut2 * (vc[2] as u64))
            >> (16 + 9);
        let va = ATAN_REDUCE[i as usize];
        let ta = f64::copysign(1.0, x) * f64::from_bits(va.0);
        let id = f64::copysign(1.0, x) * i as f64;
        h = (x - ta) / (1. + x * ta);
        a = Dekker::new(
            f64::copysign(1.0, x) * f64::from_bits(va.1) + f64::from_bits(0x3c88469898cc5170) * id,
            f64::from_bits(0x3f8921fb54442d00) * id,
        );
    }
    let h2 = h * h;
    let h4 = h2 * h2;

    let f0 = f_fmla(h2, f64::from_bits(CH[3]), f64::from_bits(CH[2]));
    let f1 = f_fmla(h2, f64::from_bits(CH[1]), f64::from_bits(CH[0]));

    let f = f_fmla(h4, f0, f1);
    a.lo = f_fmla(h, f, a.lo);
    // begin_atanpi
    /* Now ah + al approximates atan(x) with error bounded by 0x3.fp-52*h
    (see atan.c), thus by 0x1.41p-52*h after multiplication by 1/pi.
    We normalize ah+al so that the rounding error in muldd is negligible
    below. */
    let e0 = h * f64::from_bits(0x3ccf800000000000); // original value in atan.c
    let ub0 = (a.lo + e0) + a.hi; // original value in atan.c
    a = Dekker::from_exact_add(a.hi, a.lo);
    a = Dekker::quick_mult(a, Dekker::new(ONE_OVER_PIL, ONE_OVER_PIH));
    /* The rounding error in muldd() and the approximation error between 1/pi
    and ONE_OVER_PIH+ONE_OVER_PIL are absorbed when rounding up 0x3.fp-52*pi
    to 0x1.41p-52. */
    let e = h * f64::from_bits(0x3cb4100000000000); // atanpi_specific
    // end_atanpi
    let ub = (a.lo + e) + a.hi;
    let lb = (a.lo - e) + a.hi;
    if ub == lb {
        return ub;
    }
    as_atan_refine2(x, ub0)
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn atanpi_test() {
        assert_eq!(0.000014571070806516354, f_atanpi(0.00004577636903266291));
        assert_eq!(-0.000014571070806516354, f_atanpi(-0.00004577636903266291));
        assert_eq!(-0.13664770469904508, f_atanpi(-0.4577636903266291));
    }
}
