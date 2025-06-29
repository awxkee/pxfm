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
use crate::acospi_table::{ACOSPI_ERR, ACOSPI_TABLE};
use crate::common::f_fmla;
use crate::dekker::Dekker;

const PI_HI: f64 = f64::from_bits(0x400921fb54442d18);
const PI_LO: f64 = f64::from_bits(0x3ca1a62633145c07);

const ONE_OVER_PIH: f64 = f64::from_bits(0x3fd45f306dc9c883);
const ONE_OVER_PIL: f64 = f64::from_bits(0xbc76b01ec5417056);

/// Computes acos(x)/PI
///
/// Max ULP 1.5 on |x| < 0.99999
#[inline]
pub fn f_acospi(x: f64) -> f64 {
    let u = x.to_bits() & 0x7fff_ffff_ffff_ffff;
    let absx = f64::from_bits(u);
    let v_x: u64;
    let u_bytes = u.to_ne_bytes();
    let k = i32::from_ne_bytes([u_bytes[4], u_bytes[5], u_bytes[6], u_bytes[7]]);
    let u_low = u32::from_ne_bytes([u_bytes[0], u_bytes[1], u_bytes[2], u_bytes[3]]);
    if k < 0x3fe80000 {
        /* |x| < 0.75 */
        // avoid spurious underflow:
        // for |x| <= 0x1.921fb54442d18p-54, acospi(x) rounds to 0.5 to nearest
        if k < 0x3c9921fb {
            // acospi(x) ~ 1/2 - x/pi
            return f_fmla(f64::from_bits(0xbc80000000000000), x, 0.5);
        }
        /* approximate acos(x) by p(x-xmid), where [0,0.75) is split
        into 192 sub-intervals */
        v_x = (1.0 + absx).to_bits(); /* 1 <= v.x < 2 */
        /* v.i[HIGH] contains 20 significant bits in its low bits, we shift by 12
        to get the upper 8 (ignoring the implicit leading bit) */
        let mut i = ((v_x >> (12 + 32)) & 255) as i32;
        if i == 192 {
            i = 191;
        }
        let p = ACOSPI_TABLE[i as usize];
        let y = absx - f64::from_bits(p.7); /* p[7] = xmid */

        let yy = y * y;
        /* evaluate in parallel p[1] + p[2] * y and p[3] + p[4] * y, and
        p[5] + p[6] * y using Estrin's scheme */
        let p56 = f_fmla(f64::from_bits(p.6), y, f64::from_bits(p.5));
        let p34 = f_fmla(f64::from_bits(p.4), y, f64::from_bits(p.3));
        let mut zh = f_fmla(p56, yy, p34);
        zh = f_fmla(zh, y, f64::from_bits(p.2));
        let z = Dekker::from_exact_add(f64::from_bits(p.1), y * zh);
        let mut d = Dekker::from_exact_add(f64::from_bits(p.0), z.hi * y);
        d.lo = f_fmla(z.lo, y, d.lo);
        /* Special case for i=0, since we are obliged to use xmid=0 (so that
        x-xmid is exact) thus we can't use Gal's trick.  This costs about
        0.5 cycle in the average time (for both branches).  */
        if i == 0 {
            d.hi += f64::from_bits(0x3c91a62792c17e8c);
        }
        /* acos(x) ~ du + dv for x > 0, pi - (u + v) for x < 0 */
        if x < 0.
        /* acos(-x) = pi-acos(x) */
        {
            let p = Dekker::from_exact_add(PI_HI, -d.hi);
            d.hi = p.hi;
            d.lo = PI_LO + p.lo - d.lo;
        }

        // acospi_begin
        /* We multiply the approximation u+v, with maximal error say 2^-e
        by 1/pi. The maximal value of |u+v| is less than 2.42 (for x=-0.75).
        The maximal error is the sum of several terms:
        * 2^-e * (ONE_OVER_PIH + ONE_OVER_PIL) < 2^-e * 2^-1.651
        * (u+v)*|ONE_OVER_PIH+ONE_OVER_PIL-1/pi| < 2.42*2^-109.523 < 2^-108
        * the ignored term v*ONE_OVER_PIL in d_mul. The maximal observed value
          of v is 0x1.06d413839cafcp-51 for x=-0x1.6a01f2fb71p-1 (rndd),
          we conjecture |v| < 2^-50. Then |v*ONE_OVER_PIL| < 2^-105
        * the rounding error in d_mul. The d_mul call decomposes into:
          a_mul (u, lo1, u_in, ONE_OVER_PIH)
          lo2 = __builtin_fma (u_in, ONE_OVER_PIL, lo1)
          v = __builtin_fma (v_in, ONE_OVER_PIH, lo2)
          Since |u| <= acos(-0.75)/pi < 0.8 we have |lo1| <= ulp(0.8) <= 2^-53.
          Then since |u_in| <= 2.42, |lo2| <= |2.42*ONE_OVER_PIL|+2^-53
                                           < 2^-52.485
          Then |v| <= 2^-50+ONE_OVER_PIH*2^-52.485 < 2^-49.920.
          The rounding error is bounded by ulp(lo2)+ulp(v) <= 2^-105+2^-102
          < 2^-101.83.
        The total error is thus bounded by:
        2^-e * 2^-1.651 + 2^-108 + 2^-105 + 2^-101.83 < Err[i]
        */

        d = Dekker::quick_mult(d, Dekker::new(ONE_OVER_PIL, ONE_OVER_PIH));

        // acospi_end

        let err = ACOSPI_ERR[i as usize]; // acospi_specific

        d.hi + (d.lo - f64::from_bits(err))
        // right = du + (dv + err);
        // if (__builtin_expect (left != right, 0))
        // return accurate_path (x); /* hard to round case */
        // return left;
    }
    /*--------------------------- 0.75 <= |x| < 1 ---------------------*/
    else if k < 0x3ff00000 {
        /* |x| < 1 */
        /* approximate acos(x) by sqrt(1-x)*p(x-xmid) where p is a polynomial,
        and [0.75,1) is split into 64 sub-intervals */
        v_x = (1.0 + absx).to_bits(); /* 1 <= v.x <= 2 */
        /* The low 20 bits of v.i[HIGH] are the upper bits (except the
        implicit leading bit) of the significand of 1+|x|.
        Warning: v.x might be 2 for rounding up or nearest. */
        let i = if f64::from_bits(v_x) == 2.0 {
            255
        } else {
            ((v_x >> 32) & 0xff000) >> 12
        };
        let p = ACOSPI_TABLE[i as usize];
        let y = absx - f64::from_bits(p.6); /* exact (p[6] = xmid) */
        let h1 = 1.0 - absx; /* exact since |x| >= 0.5 */
        let dh1 = Dekker::from_exact_sqrt(h1);
        /* use Estrin's scheme to evaluate p2 + p3*y + p4*y^2 + p5*y^3 */
        let yy = y * y;
        let p45 = f_fmla(f64::from_bits(p.5), y, f64::from_bits(p.4));
        let p23 = f_fmla(f64::from_bits(p.3), y, f64::from_bits(p.2));
        let mut zh = f_fmla(p45, yy, p23);
        zh = f_fmla(zh, y, f64::from_bits(p.1));
        let mut z = Dekker::from_exact_add(f64::from_bits(p.0), zh * y);
        let l1zh = dh1.lo * z.hi; /* compute earlier */
        let h1zl = dh1.hi * z.lo;
        /* acos(x) ~ (h1 + l1) * (zh + zl) */
        let mut vd = Dekker::from_exact_mult(dh1.hi, z.hi);
        vd.lo += l1zh + h1zl;
        if x < 0.
        /* acos(x) = pi - (u+v) */
        {
            // let jk = Dekker::from_exact_add (&du, &zl, pi_hi, -du);
            let jk = Dekker::from_exact_add(PI_HI, -vd.hi);
            /* acos(x) = u + zl + pi_lo - v */
            vd.lo = z.lo + PI_LO - vd.lo;
            z.lo = jk.lo;
            vd.hi = jk.hi;
        }

        // acospi_begin
        /* Similar analysis as above.
        We multiply the approximation u+v, with maximal error 2^-e
        by 1/pi. The maximal value of |u+v| is pi (for x=-1).
        The maximal error is the sum of several terms:
        * 2^-e * (ONE_OVER_PIH + ONE_OVER_PIL) < 2^-e * 2^-1.651
        * (u+v)*|ONE_OVER_PIH+ONE_OVER_PIL-1/pi| < pi*2^-109.523 < 2^-107
        * the ignored term v*ONE_OVER_PIL in d_mul. The maximal observed value
          of v is 0x1.4586d502c6913p-51 for x=-0x1.fffcc87dece8p-1 (rndd),
          we conjecture |v| < 2^-50. Then |v*ONE_OVER_PIL| < 2^-105
        * the rounding error in d_mul. The d_mul call decomposes into:
          a_mul (u, lo1, u_in, ONE_OVER_PIH)
          lo2 = __builtin_fma (u_in, ONE_OVER_PIL, lo1)
          v = __builtin_fma (v_in, ONE_OVER_PIH, lo2)
          Since |u| <= acospi(-1) < 1 we have |lo1| <= ulp(1-) <= 2^-53.
          Then since |u_in| <= pi, |lo2| <= |pi*ONE_OVER_PIL|+2^-53
                                           < 2^-52.361.
          Then |v| <= 2^-50+ONE_OVER_PIH*2^-52.361 < 2^-49.913.
          The rounding error is bounded by ulp(lo2)+ulp(v) <= 2^-105+2^-102
          < 2^-101.83.
        The total error is thus bounded by:
        2^-e * 2^-1.651 + 2^-107 + 2^-105 + 2^-101.83 < Err[i]
        */
        let d = Dekker::quick_mult(vd, Dekker::new(ONE_OVER_PIL, ONE_OVER_PIH));
        // acospi_end

        let err = ACOSPI_ERR[i as usize]; // acospi_specific
        d.hi + (d.lo - f64::from_bits(err))
    }
    /*   else  if (k < 0x3ff00000)    */

    /*---------------------------- |x|>=1 -----------------------*/
    else if k == 0x3ff00000 && u_low == 0 {
        if x > 0. { 0. } else { 1. }
    }
    // acospi_specific
    else if k > 0x7ff00000 || (k == 0x7ff00000 && u_low != 0) {
        x + x // case x=nan
    } else {
        f64::NAN
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn acospi_test() {
        assert_eq!(f_acospi(0.5), 0.3333333333333333);
        assert!(f_acospi(1.5).is_nan());
    }
}
